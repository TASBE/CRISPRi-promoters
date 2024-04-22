import os
import sbol3
import latex_generation
from shared_global_names import MODEL_FILE

# Set the working directory to be the SBOL sub-folder
os.chdir('sbol')

doc = sbol3.Document()
print(f'Reading {MODEL_FILE}')
doc.read(MODEL_FILE)

# Generate a table of symbols in its own document
with open('../equations/generated_table.tex', 'w') as out:
    print('Writing table of symbols')
    out.write(latex_generation.make_symbol_table())

# For each system in the document, write the equations
with open('../equations/generated_equations.tex', 'w') as out:
    # Get all of the systems and sort by the display name (so that 0 and 1
    # target sites come first and so that the heterogeneous and identical 
    # models are interleaved)
    system_components = [o for o in doc.objects if isinstance(o, sbol3.Component)]
    system_components.sort(key=lambda x: x.display_name)
    # Loop through the systems and write the equations
    for c in system_components:
        print(f'Writing model for {c.identity}')
        out.write(latex_generation.make_latex_model(c))
