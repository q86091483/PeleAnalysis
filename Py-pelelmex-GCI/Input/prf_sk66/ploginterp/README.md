# PLOG reaction pressure interpolation

This python utility interpolates all the PLOG reactions to a constant pressure in a mechanism. This is useful to later convert the mechanism to the format used in the Pele combustion suite. 

## Usage

1. Clone the repository in the file where the Cantera ‘.yaml‘ mechanism is

‘git clone https://git.step.polymtl.ca/tamarin/ploginterp.git‘

2. Install the requirements (preferably in a virtual environment)

‘‘‘bash
cd ploginterp
pip install -f requirements.txt
‘‘‘

3. Run the conversion utility on your mechanism

For example convert my_mechanism.yaml to a constant 2.0 atm pressure

‘‘‘
cd ..
python -m ploginterp my_mechanism.yaml 2.0
‘‘‘
The file my_mechanism_10_atm.yaml should be created 
