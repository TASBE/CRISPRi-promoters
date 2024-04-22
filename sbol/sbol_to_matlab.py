import os

import sbol3

import matlab_generation
from shared_global_names import MODEL_FILE

# Set the working directory to be the SBOL sub-folder
# Needed when starting from the run.debug button in the IDE
os.chdir('sbol')


doc = sbol3.Document()
print(f'Reading {MODEL_FILE}')
doc.read(MODEL_FILE)

# For each system in the document, generate a matlab model
# Model has three parts:
for c in (o for o in doc.objects if isinstance(o, sbol3.Component)):
    with open(os.path.join('../models', f'{c.display_id}.m'),
              'w') as out:
        print(f'Writing model for {c.identity}')
        model, parameters = matlab_generation.make_matlab_model(c, 'ode15s')
        out.write(model)
