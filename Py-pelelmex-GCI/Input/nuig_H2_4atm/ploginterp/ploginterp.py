import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy.optimize import curve_fit


def extract_arrhenius_params(rate):
    """ Extract the parameters from a Cantera Arrhenius rate """
    A = rate.pre_exponential_factor
    b = rate.temperature_exponent
    E = rate.activation_energy
    return A, b, E


def Arrhenius_evaluate(T, A, b, E):
    """ Compute the Arrhenius reaction rate from a temperature value and parameters """
    k = A * T ** b * np.exp(E / (ct.gas_constant * T))
    return k


def pressure_interpolate_reaction(reaction, pressure):
    """
    Interpolate the Arrhenius rate for a specific pressure (P)
    :param pressure: Fixed pressure at wich the mechanism is interpolated (atm)
    """
    # For convenience 
    log = np.log
    # Conversion to Pa
    P = pressure * ct.one_atm
    # Extract the pressure array
    pressures = np.array([ri[0] for ri in reaction.rate.rates])
    # Extract the Arrhenius rate array
    rates = np.array([ri[1] for ri in reaction.rate.rates])
    
    # Sort the pressure if needed
    if ~np.all(pressures[:-1] <= pressures[1:]):
        rates = rates[np.argsort(pressures)]
        pressures = np.sort(pressures)
        
    # Check if the interpolation pressure is below the lowest pressure
    if P < pressures[0]:
        A, b, E = extract_arrhenius_params(rates[0])

    # Check if the interpolation pressure is above the highest pressure
    elif P > pressures[-1]:
        A, b, E = extract_arrhenius_params(rates[-1])

    # The pressure is between two available presssures points
    else:
        
        # Check if the pressure is less than 1% different from a point
        P_diff = np.abs(pressures - P)
        
        if np.min(P_diff) < (pressures * 0.01)[np.argmin(P_diff)]:
            A, b, E = extract_arrhenius_params(rates[np.argmin(P_diff)])
        
        # If not we interpolate
        else:
        
            # Find the index of the two pressures to be used in the interpolation
            idx1 = np.arange(pressures.size)[pressures <= P][-1] 
            idx2 = np.arange(pressures.size)[pressures >= P][0]
        
            # Get the two pressure values
            P1 = pressures[idx1]
            P2 = pressures[idx2]

            # Get the arrhenius rate values
            A1, b1, E1 = extract_arrhenius_params(rates[idx1])
            A2, b2, E2 = extract_arrhenius_params(rates[idx2])
            
            # Scale the b exponent values so they are positive
            min_b = np.min([b1, b2])
            b_scaled = np.array([b1, b2]) + 1 - min_b
            # Compute guess values with a log interpolation
            A_guess = np.exp(np.interp(log(P), log([P1, P2]), log([A1, A2])))
            b_guess = np.exp(np.interp(log(P), log([P1, P2]), log([b_scaled[0], b_scaled[1]])))
            E_guess = np.exp(np.interp(log(P), log([P1, P2]), log([E1, E2])))
            # Scale back the b values
            b_guess = b_guess - 1 + min_b

            # Create a temperature vector for the curve fitting
            T = np.linspace(300, 2500, 100)
            
            # Compute the reaction rates at the two pressure points
            k1 = Arrhenius_evaluate(T, A1, b1, E1)
            k2 = Arrhenius_evaluate(T, A2, b2, E2)
            
            # Compute the interpolated reaction rate value
            # From the Cantera documentation : 
            # https://cantera.org/documentation/docs-2.4/doxygen/html/dd/d89/classCantera_1_1Plog.html
            k = np.exp(log(k1) + (log(k2) - log(k1)) * ((log(P) - log(P1)) / (log(P2) - log(P1))))
            
            # Fit the A, b, E parameters to the interpolated curve 
            (A, b, E), _ = curve_fit(Arrhenius_evaluate, T, k, p0=[A_guess, b_guess, E_guess])
        
    return A, b, E


def validate_temp_point(plog_reactions, conv_reactions, temperature, pressure):
    """
    Validate the reaction rates between plog reactions and interpolated reactions
    """
    print(f'\n Variable pressure and interpolated fixed pressure reaction rates at {temperature}K \n')
    for r1, r2 in zip(plog_reactions, conv_reactions):
        logr = r1.rate(temperature, pressure * ct.one_atm)
        fixr = r2.rate(temperature)
        change = np.abs(fixr - logr) / logr
        print(f'Plog rate : {logr:<20} \t Fixed P rate : {fixr:<20} \t Relative change :  {change:<20}')    
        
        
def reformat_duplicates(mech_file, gas):
    """
    Reformat the .yaml file so that duplicate reactions are taken care of
    """
    # Find out which reactions are duplicate
    equations = [ri.equation for ri in gas.reactions()]
    is_duplicate = np.array([False for ri in gas.reactions()])

    for eq_name in equations:
        name_comp = np.where([eq_name == eqi for eqi in equations])
        if name_comp[0].size > 1:
            is_duplicate[name_comp[0]] = True

    duplicate_equations = np.array(equations)[is_duplicate]
    
    file = open('temp.yaml')
    new_file = ""
    for i, line in enumerate(file):
        new_file += line
        if np.any([di in line for di in duplicate_equations]):
            new_file += '    duplicate: true\n'
    file.close()
    fout = open(mech_file, "w")
    fout.write(new_file)
    fout.close()
    os.remove('temp.yaml')
    print(f'\n New mechanism file saved at {mech_file}')
    
    
def write_yaml_with_duplicates(input_file, gas, pressure):
    """
    Write the Solution instace to a yaml file with duplicates
    """
    gas.write_yaml('temp.yaml', skip_user_defined=True)
    pressure_string = '_' + ''.join(str(pressure).split('.')) + 'atm'
    out_file = ".".join([input_file.split(".")[0] + pressure_string, 'yaml'])
    reformat_duplicates(out_file, gas)
    

def convert_mech_to_fixed_pressure(input_file, pressure, validate=False, graphic_validate=False):
    """
    Convert a cantera kinetics mechanism with varied pressure reactions 
    to a fixed pressure value
    :param input_file: mechanism file (str)
    :param pressure: pressure used for the interpolation in atm (float) 
    :param validate: Perform a validation of the interpolated reactions (bool)
    :param graphic_validate: Perform a graphic validation 
    of the interpolated reactions (bool)
    
    :return : Cantera Solution instance of the new mechanism
    """
    # Load the reactions from the original file
    ref_gas = ct.Solution(input_file)
    species = ct.Species.list_from_file(input_file)
    reference_phase = ct.Solution(thermo='ideal-gas', kinetics='gas', species=species)
    all_reactions = ct.Reaction.list_from_file(input_file, reference_phase)
    # Arrays for the converted reactions
    new_reactions = []
    plog_reactions = []
    conv_reactions = []
    # Rate types wich do not need to be converted
    legal_rates_types = ['Arrhenius', 'Troe', 'Lindemann', 'falloff']
    for R in all_reactions:
        # If not pressure dependent
        if R.rate.type in legal_rates_types:
            # Append to the new reactions
            new_reactions.append(R)
        # Else convert to arrhenius
        else:
            print(R.rate.type)
            # Original reactions that were converted
            plog_reactions.append(R)
            # Fixed pressure Arrhenius values 
            A, b, E = pressure_interpolate_reaction(R, pressure)
            # Instantiate a new Cantera reaction object
            new_R = ct.Reaction(equation=R.equation, rate=ct.ArrheniusRate(A, b, E))
            # Append to the new reactions and store created reactions objects for validation
            new_reactions.append(new_R)
            conv_reactions.append(new_R)
            
    new_gas = ct.Solution(name="fixed_pressure_mech", thermo=ref_gas.thermo_model, kinetics="gas", 
                          transport_model=ref_gas.transport_model, species=species, reactions=new_reactions)
    
    if validate:
        print('Reactions converted to fixed pressure : \n')
        for ri in plog_reactions:
            print(ri.equation)
        validate_temp_point(plog_reactions, conv_reactions, 300, pressure)
        validate_temp_point(plog_reactions, conv_reactions, 1000, pressure)
        validate_temp_point(plog_reactions, conv_reactions, 3000, pressure)
        
    if graphic_validate:
        # Graphical validation
        T = np.linspace(300, 3000, 200)
        graph_len = len(plog_reactions) // 3
        extra_len = np.min([len(plog_reactions) % 3, 1])
        fig, axs = plt.subplots(graph_len + extra_len, 3, figsize=(14, 30))

        for r1, r2, ax in zip(plog_reactions, conv_reactions, axs.flatten()):
            plt.sca(ax)
            k1 = np.array([r1.rate(ti, 0.7 * ct.one_atm) for ti in T])
            k2 = np.array([r2.rate(ti) for ti in T])
            plt.plot(T, k1, color='r', label='Plog', alpha=0.9)
            plt.plot(T, k2, color='b', linestyle='--', label='Interp', alpha=0.9)
            plt.title(f'{r1.equation}')
            plt.legend()
            plt.grid()
            ax.set_ylabel('Reaction rate')
            ax.set_xlabel('Temperature (K)')
        plt.tight_layout()
        fig.savefig('Graphic_validation', dpi=500)
    
    write_yaml_with_duplicates(input_file, new_gas, pressure)
    

file = sys.argv[1]
pressure = float(sys.argv[2])

convert_mech_to_fixed_pressure(file, pressure, validate=True)
