# Import the necessary packages
import sbol3
from sbol_utilities.component import add_feature, add_interaction

# Import functions that make the circuit modules
import builders

# Set global variables
MODEL_FILE = 'sbol/gRNA_models.nt'  # Assumes running from the root directory
PROJECT_NAMESPACE = 'http://bbn.com/apt-dcas9-regulation'

# Make the document, set the namespace
doc = sbol3.Document()
sbol3.set_namespace(PROJECT_NAMESPACE)

###############################################################################
# Constitutive Expression
###############################################################################
system = sbol3.Component('No_gRNA_control', sbol3.SBO_FUNCTIONAL_ENTITY,
                         name="0 Target Sites")
doc.add(system)
# Create a GFP-CRISPRa system
vector_2 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                                                       name='V2'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_2: sbol3.SBO_REACTANT})
dCas9_pro, gfp_promoter, gfp_pro = builders.make_gfp_dCas9_module(vector_2)
# Set the GFP expression as constitutive
builders.constitutive(gfp_promoter)
# Set the inputs and outputs as an interface, for generating the ODEs
# TODO: Warning will go away after resolution of https://github.com/SynBioDex/pySBOL3/issues/324
inputs = [vector_2]
system.interface = sbol3.Interface(inputs=inputs, outputs=[gfp_pro])

###############################################################################
# 1 Target Site
###############################################################################
system = sbol3.Component('Single_gRNA_repression',
                         sbol3.SBO_FUNCTIONAL_ENTITY,
                         name="1 Target Site")
doc.add(system)
# Create a GFP-CRISPRa system
vector_2 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                       name='V2'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_2: sbol3.SBO_REACTANT})
dCas9_pro, gfp_promoter, gfp_pro = builders.make_gfp_dCas9_module(vector_2)
# Create a gRNA system on a separate vector
vector_1 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                       name='V1'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_1: sbol3.SBO_REACTANT})
dCas9_gRNA1 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        1)
# Set the inputs and outputs as an interface, for generating the ODEs
# TODO: Warning will go away after resolution of https://github.com/SynBioDex/pySBOL3/issues/324
inputs = [vector_1, vector_2]
system.interface = sbol3.Interface(inputs=inputs, outputs=[gfp_pro])

###############################################################################
# 2 Heterogeneous Target Sites
###############################################################################
system = sbol3.Component('Multiplexed_2_gRNA_Repression',
                         sbol3.SBO_FUNCTIONAL_ENTITY,
                         name="2 Heterogeneous Target Sites")
doc.add(system)
# Create a GFP-CRISPRa system
vector_2 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                                                       name='V2'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_2: sbol3.SBO_REACTANT})
dCas9_pro, gfp_promoter, gfp_pro = builders.make_gfp_dCas9_module(vector_2)
# Create two gRNA systems on a separate vector
vector_1 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                                                       name='V1'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_1: sbol3.SBO_REACTANT})
dCas9_gRNA1 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=1)
dCas9_gRNA2 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=2)
# Assuming there is no interaction between the gRNAs, if there were use the
# builders.add_complex_interference function to add that
# Set the inputs and outputs as an interface, for generating the ODEs
# TODO: Warning will go away after resolution of https://github.com/SynBioDex/pySBOL3/issues/324
inputs = [vector_1, vector_2]
system.interface = sbol3.Interface(inputs=inputs, outputs=[gfp_pro])

###############################################################################
# 3 Heterogeneous Target Sites
###############################################################################
system = sbol3.Component('Multiplexed_3_gRNA_Repression',
                         sbol3.SBO_FUNCTIONAL_ENTITY,
                         name="3 Heterogeneous Target Sites")
doc.add(system)
# Create a GFP-CRISPRa system
vector_2 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                                                       name='V2'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_2: sbol3.SBO_REACTANT})
dCas9_pro, gfp_promoter, gfp_pro = builders.make_gfp_dCas9_module(vector_2)
# Create two gRNA systems on a separate vector
vector_1 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                                                       name='V1'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_1: sbol3.SBO_REACTANT})
dCas9_gRNA1 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=1)
dCas9_gRNA2 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=2)
dCas9_gRNA3 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=3)
# Assuming there is no interaction between the gRNAs, if there were use the
# builders.add_complex_interference function to add that
# Set the inputs and outputs as an interface, for generating the ODEs
# TODO: Warning will go away after resolution of https://github.com/SynBioDex/pySBOL3/issues/324
inputs = [vector_1, vector_2]
system.interface = sbol3.Interface(inputs=inputs, outputs=[gfp_pro])

###############################################################################
# 4 Heterogeneous Target Sites
###############################################################################
system = sbol3.Component('Multiplexed_4_gRNA_Repression',
                         sbol3.SBO_FUNCTIONAL_ENTITY,
                         name="4 Heterogeneous Target Sites")
doc.add(system)
# Create a GFP-CRISPRa system
vector_2 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                                                       name='V2'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_2: sbol3.SBO_REACTANT})
dCas9_pro, gfp_promoter, gfp_pro = builders.make_gfp_dCas9_module(vector_2)
# Create two gRNA systems on a separate vector
vector_1 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                                                       name='V1'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_1: sbol3.SBO_REACTANT})
dCas9_gRNA1 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=1)
dCas9_gRNA2 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=2)
dCas9_gRNA3 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=3)
dCas9_gRNA4 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=4)
# Assuming there is no interaction between the gRNAs, if there were use the
# builders.add_complex_interference function to add that
# Set the inputs and outputs as an interface, for generating the ODEs
# TODO: Warning will go away after resolution of https://github.com/SynBioDex/pySBOL3/issues/324
inputs = [vector_1, vector_2]
system.interface = sbol3.Interface(inputs=inputs, outputs=[gfp_pro])

###############################################################################
# 5 Heterogeneous Target Sites
###############################################################################
system = sbol3.Component('Multiplexed_5_gRNA_Repression',
                         sbol3.SBO_FUNCTIONAL_ENTITY,
                         name="5 Heterogeneous Target Sites")
doc.add(system)
# Create a GFP-CRISPRa system
vector_2 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                                                       name='V2'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_2: sbol3.SBO_REACTANT})
dCas9_pro, gfp_promoter, gfp_pro = builders.make_gfp_dCas9_module(vector_2)
# Create two gRNA systems on a separate vector
vector_1 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                                                       name='V1'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_1: sbol3.SBO_REACTANT})
dCas9_gRNA1 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=1)
dCas9_gRNA2 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=2)
dCas9_gRNA3 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=3)
dCas9_gRNA4 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=4)
dCas9_gRNA5 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=5)
# Assuming there is no interaction between the gRNAs, if there were use the
# builders.add_complex_interference function to add that
# Set the inputs and outputs as an interface, for generating the ODEs
# TODO: Warning will go away after resolution of https://github.com/SynBioDex/pySBOL3/issues/324
inputs = [vector_1, vector_2]
system.interface = sbol3.Interface(inputs=inputs, outputs=[gfp_pro])

###############################################################################
# 6 Heterogeneous Target Sites
###############################################################################
system = sbol3.Component('Multiplexed_6_gRNA_Repression',
                         sbol3.SBO_FUNCTIONAL_ENTITY,
                         name="6 Heterogeneous Target Sites")
doc.add(system)
# Create a GFP-CRISPRa system
vector_2 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                                                       name='V2'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_2: sbol3.SBO_REACTANT})
dCas9_pro, gfp_promoter, gfp_pro = builders.make_gfp_dCas9_module(vector_2)
# Create two gRNA systems on a separate vector
vector_1 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                                                       name='V1'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_1: sbol3.SBO_REACTANT})
dCas9_gRNA1 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=1)
dCas9_gRNA2 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=2)
dCas9_gRNA3 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=3)
dCas9_gRNA4 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=4)
dCas9_gRNA5 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=5)
dCas9_gRNA6 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=6)
# Assuming there is no interaction between the gRNAs, if there were use the
# builders.add_complex_interference function to add that
# Set the inputs and outputs as an interface, for generating the ODEs
# TODO: Warning will go away after resolution of https://github.com/SynBioDex/pySBOL3/issues/324
inputs = [vector_1, vector_2]
system.interface = sbol3.Interface(inputs=inputs, outputs=[gfp_pro])

###############################################################################
# 2 Identical Target Sites
###############################################################################
system = sbol3.Component('Multisite_2_gRNA_Repression',
                         sbol3.SBO_FUNCTIONAL_ENTITY,
                         name="2 Identical Target Sites")
doc.add(system)

# Create a GFP-CRISPRa system
vector_2 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                                                       name='V2'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_2: sbol3.SBO_REACTANT})
dCas9_pro, gfp_promoter, gfp_pro = builders.make_gfp_dCas9_module(vector_2)

# Create one gRNA systems on a separate vector
vector_1 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                                                       name='V1'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_1: sbol3.SBO_REACTANT})
dCas9_gRNA1 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=1,
                                        n_interactions=2)

# Set the inputs and outputs as an interface, for generating the ODEs
# TODO: Warning will go away after resolution of https://github.com/SynBioDex/pySBOL3/issues/324
inputs = [vector_1, vector_2]
system.interface = sbol3.Interface(inputs=inputs, outputs=[gfp_pro])

###############################################################################
# 3 Identical Target Sites
###############################################################################
system = sbol3.Component('Multisite_3_gRNA_Repression',
                         sbol3.SBO_FUNCTIONAL_ENTITY,
                         name="3 Identical Target Sites")
doc.add(system)

# Create a GFP-CRISPRa system
vector_2 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                                                       name='V2'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_2: sbol3.SBO_REACTANT})
dCas9_pro, gfp_promoter, gfp_pro = builders.make_gfp_dCas9_module(vector_2)

# Create one gRNA systems on a separate vector
vector_1 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                                                       name='V1'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_1: sbol3.SBO_REACTANT})
dCas9_gRNA1 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=1,
                                        n_interactions=3)

# Set the inputs and outputs as an interface, for generating the ODEs
# TODO: Warning will go away after resolution of https://github.com/SynBioDex/pySBOL3/issues/324
inputs = [vector_1, vector_2]
system.interface = sbol3.Interface(inputs=inputs, outputs=[gfp_pro])

###############################################################################
# 4 Identical Target Sites
###############################################################################
system = sbol3.Component('Multisite_4_gRNA_Repression',
                         sbol3.SBO_FUNCTIONAL_ENTITY,
                         name="4 Identical Target Sites")
doc.add(system)

# Create a GFP-CRISPRa system
vector_2 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                                                       name='V2'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_2: sbol3.SBO_REACTANT})
dCas9_pro, gfp_promoter, gfp_pro = builders.make_gfp_dCas9_module(vector_2)

# Create one gRNA systems on a separate vector
vector_1 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                                                       name='V1'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_1: sbol3.SBO_REACTANT})
dCas9_gRNA1 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=1,
                                        n_interactions=4)

# Set the inputs and outputs as an interface, for generating the ODEs
# TODO: Warning will go away after resolution of https://github.com/SynBioDex/pySBOL3/issues/324
inputs = [vector_1, vector_2]
system.interface = sbol3.Interface(inputs=inputs, outputs=[gfp_pro])

###############################################################################
# 5 Identical Target Sites
###############################################################################
system = sbol3.Component('Multisite_5_gRNA_Repression',
                         sbol3.SBO_FUNCTIONAL_ENTITY,
                         name="5 Identical Target Sites")
doc.add(system)

# Create a GFP-CRISPRa system
vector_2 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                                                       name='V2'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_2: sbol3.SBO_REACTANT})
dCas9_pro, gfp_promoter, gfp_pro = builders.make_gfp_dCas9_module(vector_2)

# Create one gRNA systems on a separate vector
vector_1 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                                                       name='V1'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_1: sbol3.SBO_REACTANT})
dCas9_gRNA1 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=1,
                                        n_interactions=5)

# Set the inputs and outputs as an interface, for generating the ODEs
# TODO: Warning will go away after resolution of https://github.com/SynBioDex/pySBOL3/issues/324
inputs = [vector_1, vector_2]
system.interface = sbol3.Interface(inputs=inputs, outputs=[gfp_pro])

###############################################################################
# 6 Identical Target Sites
###############################################################################
system = sbol3.Component('Multisite_6_gRNA_Repression',
                         sbol3.SBO_FUNCTIONAL_ENTITY,
                         name="6 Identical Target Sites")
doc.add(system)

# Create a GFP-CRISPRa system
vector_2 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                                                       name='V2'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_2: sbol3.SBO_REACTANT})
dCas9_pro, gfp_promoter, gfp_pro = builders.make_gfp_dCas9_module(vector_2)

# Create one gRNA systems on a separate vector
vector_1 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA],
                                                       name='V1'))
add_interaction(sbol3.SBO_DEGRADATION,
                name='Vector degradation',
                participants={vector_1: sbol3.SBO_REACTANT})
dCas9_gRNA1 = builders.make_gRNA_module(vector_1,
                                        dCas9_pro,
                                        gfp_promoter,
                                        True,
                                        idx_gRNA=1,
                                        n_interactions=6)

# Set the inputs and outputs as an interface, for generating the ODEs
# TODO: Warning will go away after resolution of https://github.com/SynBioDex/pySBOL3/issues/324
inputs = [vector_1, vector_2]
system.interface = sbol3.Interface(inputs=inputs, outputs=[gfp_pro])

###############################################################################
# Write the model file #
###############################################################################
doc.write(MODEL_FILE, sbol3.SORTED_NTRIPLES)
