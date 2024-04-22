import logging
import itertools
from typing import Dict, List, Optional
import re

import sbol3
import tyto
from sbol_utilities.helper_functions import id_sort
from sbol_utilities.component import all_in_role, in_role

from shared_global_names import RECOMBINATION

# TODO: Sort this dictionary however I want it to be in the table
# Alphabetical? Latin first then greek?
name_to_symbol = {
    'V1': ('\\constructGen{1}',
           'Construct vector, species 1'),
    'V2': ('\\constructGen{2}',
           'Construct vector, species 1'),
    'gRNA1': ('\\gRna{1}',
              'gRNA, species 1'),
    'gRNA2': ('\\gRna{2}',
              'gRNA, species 2'),
    'gRNA3': ('\\gRna{3}',
              'gRNA, species 3'),
    'gRNA4': ('\\gRna{4}',
              'gRNA, species 4'),
    'gRNA5': ('\\gRna{5}',
              'gRNA, species 5'),
    'gRNA6': ('\\gRna{6}',
              'gRNA, species 6'),
    'dCas9': ('\\proSp{\cas{}}',
              'dCas9 protein species'),
    'GFP': ('\\proSp{\gfp}',
            'GFP protein species'),
    'dCas9-gRNA1': ('\\cplx{\\cas}{1}',
                    'Complex of dCas9 protein and gRNA1'),
    'dCas9-gRNA2': ('\\cplx{\\cas}{2}',
                    'Complex of dCas9 protein and gRNA2'),
    'dCas9-gRNA3': ('\\cplx{\\cas}{3}',
                    'Complex of dCas9 protein and gRNA3'),
    'dCas9-gRNA4': ('\\cplx{\\cas}{4}',
                    'Complex of dCas9 protein and gRNA4'),
    'dCas9-gRNA5': ('\\cplx{\\cas}{5}',
                    'Complex of dCas9 protein and gRNA5'),
    'dCas9-gRNA6': ('\\cplx{\\cas}{6}',
                    'Complex of dCas9 protein and gRNA6'),
    'Transcription rate': ('\\txRate{i}',
                           'Transcription rate, species i'),
    'Transcription-translation rate': ('\\txtlRate{i}',
                                       ('Coupled transcription and translation'
                                        'rate, species i')),
    'gRNA degradation': ('\\rnaDegradeRate{}',
                         'gRNA degradation rate'),
    'Protein degradation': ('\\proDegradeRate{}',
                            'Protein degradation rate'),
    'Stable Molecule Dilution': ('\\dilutionRate{}',
                                 'Stable Molecule Dilution Rate'),
    'dCas9 degradation': ('\\proDegradeRate{\proSp{\cas{}}}',
                          'dCas9 and dCas9-gRNA complex degradation rate'),
    'dCas9-gRNA1-complex-degradation': ('\\casCompDegradeRate{}',
                                        'dCas9-gRNA1 complex degradation rate'),
    'dCas9-gRNA2-complex-degradation': ('\\casCompDegradeRate{}',
                                        'dCas9-gRNA2 complex degradation rate'),
    'dCas9-gRNA3-complex-degradation': ('\\casCompDegradeRate{}',
                                        'dCas9-gRNA3 complex degradation rate'),
    'dCas9-gRNA4-complex-degradation': ('\\casCompDegradeRate{}',
                                        'dCas9-gRNA4 complex degradation rate'),
    'dCas9-gRNA5-complex-degradation': ('\\casCompDegradeRate{}',
                                        'dCas9-gRNA5 complex degradation rate'),
    'dCas9-gRNA6-complex-degradation': ('\\casCompDegradeRate{}',
                                        'dCas9-gRNA6 complex degradation rate'),
    'GFP degradation': ('\\proDegradeRate{\proSp{\gfp}}',
                        'GFP protein degradation rate'),
    'Cas-gRNA binding': ('\\gRnaBind{}',
                         'dCas9 and gRNA binding rate'),
    'Hill Coefficient': ('\\hillCoeff', 'Hill coefficient'),
    'Hill Activation Constant': ('\\hillActivationConst',
                                 'Hill Activation Constant'),
    'Hill Repression Constant': ('\\hillRepressionConst',
                                 'Hill Repression Constant'),
    'Interference Matrix': ('\\intMatrix{i}{j}',
                            'Interference Matrix entry at coordinates i, j')
}
"""Dictionary mapping from SBOL names to LaTeX symbols in our convention"""


def make_symbol_table():
    """Scan the dictionary that has the symbols and their meanings and generate
    a table defining each symbol

    Arguments:
       None

    Returns:
        latex (string): String serialization of LaTeX table
    """
    # Make header
    latex = """\\begin{center}
    \\begin{tabular}{|c|c|}
    \\hline
    Variable & Meaning \\\\
    \\hline
    """
    # Add an entry for every name in the name_to_symbol dictionary
    for symbol_def_pair in name_to_symbol.values():
        latex += f"""
        ${symbol_def_pair[0]}$ & {symbol_def_pair[1]} \\\\"""
    # End table
    latex += """\hline
    \\end{tabular}
    \\end{center}"""

    return latex


def transitive_closure(d: dict) -> dict:
    """Interpreting a dictionary as an acyclic directed graph, create a transitive closure of all k->v edges
    For example {1:[2,5], 2:[3,5], 3:[], 4:[5], 5:[6], 6:[]}
    returns {1:[2,3,5,6], 2:[3,5,6], 3:[], 4:[5,6], 5:[6], 6:[]}

    :param d: dictionary to close
    :return: closure dictionary
    """
    closure = d
    pending = closure
    while pending:
        # find out which are leaves
        resolvable = {k for k, v in pending.items() if not v}
        if not resolvable:
            raise ValueError(f'Cannot compute closure on cycle graph {d}')
        # union all the targets one step away
        for r in resolvable:
            closure[r] = id_sort(set(closure[r]) | set(itertools.chain(*(closure[x] for x in closure[r]))))
        # remove everything that's been resolved
        pending = {k: [x for x in v if x not in resolvable] for k, v in pending.items() if k not in resolvable }
    return closure


def maybe_concentration(feature: sbol3.Feature) -> str:
    """Determine whether we are working with a concentration or a count based on type

    :param feature: Feature to be evaluated
    :return: symbol, possibly wrapped in a concentration
    """
    if not feature.name:  # if there is no name, then it's a pass-through and we don't add a symbol
        return ''
    if feature.name not in name_to_symbol:
        raise ValueError(f'No symbol known for name: "{feature.name}"')
    symbol = name_to_symbol[feature.name][0]
    if sbol3.SBO_DNA in feature.types:
        return symbol
    else:
        return f'\\conc{{{symbol}}}'


def effective_concentration(feature: sbol3.Feature, all_interactions) -> str:
    """Create a string of the effective concentration of a regulator, based on
    its concentration and the concentration of the other regulators that
    interfere with its action

    :param feature: Feature for the regulator of interest
    :param all_interactions:
    :return: symbol, possibly modulated by interferers
    """
    # FIXME: Only add the parentheses if there are some interfering regulators
    # TODO: See if I can regenerate the all_interactions dict here, rather than needing to pass it down
    symbol = '(' + maybe_concentration(feature)
    # Look for interference interactions with other regulators
    interferers = find_interfering_regulators(feature, all_interactions)
    # For every interference interaction
    for interfering_regulator in interferers:
        # Add on the interference term to the concentration symbol
        symbol += f' + \\intMatrix{{{name_to_symbol[feature.name][0]}}}{{{name_to_symbol[interfering_regulator.name][0]}}}' + maybe_concentration(interfering_regulator)
    # Close the parentheses
    symbol += ')'

    return symbol


def differential(feature: sbol3.Feature) -> str:
    """Return the "dX/dt" term of the ODE

    :param feature: feature to get a differential for
    :return: LaTeX string
    """
    return f'\\diff{{{maybe_concentration(feature)}}}{{t}}'


def regulation_term(interaction: sbol3.Interaction,
                    all_interactions) -> str:
    """Generate a term for regulation by transcription factor or recombinase

    :param interaction: Regulation interaction to serialize
    :return: LaTeX serialization
    """
    # Need i_type to see what type of regulation is happening
    i_type = interaction.types[0]
    # Get the sbol feature of the regulator
    if i_type == sbol3.SBO_INHIBITION:
        regulator = in_role(interaction, sbol3.SBO_INHIBITOR)
    elif i_type == sbol3.SBO_STIMULATION:
        regulator = in_role(interaction, sbol3.SBO_STIMULATOR)
    # Get the effective concentration of the regulator
    # If you are assuming that there is no interference between the gRNAs, the
    # effective concentration is just the concentration of the species
    effective_regulator_concentration = maybe_concentration(regulator)
    # If there is interference between the gRNAs, include the interference
    # coefficients in the effective concentration with the following line
    # effective_regulator_concentration = effective_concentration(regulator, all_interactions)
    # Make TF Equations
    if i_type == sbol3.SBO_INHIBITION:
        # TODO: Replace K and n with variables
        return f'\\frac{{(K_R)^n}}{{(K_R)^n + {effective_regulator_concentration}^n}}'
    elif i_type == sbol3.SBO_STIMULATION:
        # TODO: Replace K and n with variables
        return f'\\frac{{{maybe_concentration(regulator)}^n}}{{(K_A)^n + {maybe_concentration(regulator)}^n}}'
    # Make Cre equations
    elif i_type == RECOMBINATION:
        target = in_role(interaction, sbol3.SBO_MODIFIED)
        if any(tyto.SO.promoter.is_ancestor_of(r) for r in target.roles):  # Cre-off
            original = in_role(interaction, sbol3.SBO_REACTANT)
            return f'\\frac{{{maybe_concentration(original)}}}{{\\vectorGen{{}}}}' # TODO: replace vectorGen w. variable
        elif any(tyto.SO.terminator.is_ancestor_of(r) for r in target.roles):  # Cre-on
            recombined = in_role(interaction, sbol3.SBO_PRODUCT)
            return f'\\frac{{{maybe_concentration(recombined)}}}{{\\vectorGen{{}}}}' # TODO: replace vectorGen w. variable
        else:
            raise ValueError(f'Cannot give term for recombination on roles {target.roles} in {interaction.identity}')
    else:
        logging.warning(f'Cannot serialize regulation {interaction.identity} of type {tyto.SBO.get_term_by_uri(i_type)}')
        return ''


def find_interfering_regulators(feature: sbol3.feature, all_interactions) -> list:
    """Find all other regulators that interfere with a given regulator

    :param feature: Feature for the regulator of interest
    :param all_interactions: (TODO: SEE IF THERE IS A WAY I CAN GET RID OF THIS)
    :return: a list of sbol features that interfere with the given regulator
    """
    # Get the interactions just for this feature
    interactions = all_interactions[feature]
    # Filter the interactions to get just the "Control" interactions
    control_interactions = [interaction
                            for interaction in interactions
                            if sbol3.SBO_CONTROL in interaction.types]
    # Search those interactions for the interactions where this feature is the
    # thing being controlled (modified), and then save the modifier to a list
    interfering_regulators = []
    for interaction in control_interactions:
        for participant in interaction.participations:
            # The participant that is not our main feature and does the
            # modifying is our interferer
            if (participant.participant != feature.identity) and (sbol3.SBO_MODIFIER in participant.roles):
                # Add the feature for that interferer to our list
                # TODO: I didn't know you could access the document from a
                # lower level object like this, see if I can use that anywhere
                # to simplify the variables I have to pass down!
                interfering_regulators.append(participant.document.find(participant.participant))

    return interfering_regulators


def interaction_to_term(feature: sbol3.Feature, interaction: sbol3.Interaction,
                        all_interactions: Dict[sbol3.Feature, List[sbol3.Interaction]],
                        regulation: Dict[sbol3.Feature, List[sbol3.Interaction]],
                        containers: Dict[sbol3.Feature, List[sbol3.Feature]]) -> Optional[str]:
    """Generate an equation term for a given interaction, with respect to the included feature

    :param feature: Target of the term
    :param interaction: Interaction to get an equation for
    :param regulation: Dictionary of regulation interactions in the system
    :param containers: Dictionary of container relationships in system
    :return: LaTeX equation term
    """
    if len(interaction.types) != 1:
        raise ValueError(f'Expected 1 interaction type but found {len(interaction.types)} in {interaction.identity}')
    if len(feature.types) != 1:
        raise ValueError(f'Expected 1 feature type but found {len(feature.types)} in {feature.identity}')
    # find the participation for this feature and its role therein
    feature_participation = [p for p in interaction.participations if p.participant == feature.identity]
    if len(feature_participation) != 1:
        raise ValueError(f'Expected feature in 1 participant, but found {len(feature_participation)} in {interaction.identity}')
    if len(feature_participation[0].roles) != 1:
        raise ValueError(f'Do not know how to serialize multi-role participation {feature_participation[0]}')
    i_type = interaction.types[0]
    f_type = feature.types[0]
    role = feature_participation[0].roles[0]
    # serialized based on interaction type and role
    if i_type == sbol3.SBO_GENETIC_PRODUCTION:
        if role == sbol3.SBO_TEMPLATE:
            return None  # templates don't get equations - they are taken as regulator for products
        elif role == sbol3.SBO_PRODUCT:
            # Using the close matches function because the species names will
            # be V1 or V2 and the entry in the name_to_symbol dictionary is V
            # species = difflib.get_close_matches('V2', list(name_to_symbol.keys()))[0]
            # But then how would I keep track of the 1 or the 2?
            species = name_to_symbol[feature.name][0]
            template = in_role(interaction, sbol3.SBO_TEMPLATE)
            # Modulation is the regulation of either the template or the 
            # product
            modulators = id_sort(regulation[feature] + regulation[template])
            modulation = ''.join(regulation_term(r, all_interactions) for r in modulators)
            # context is the constraints of the template
            context = ''.join(maybe_concentration(ct) for ct in containers[template])
            if f_type == sbol3.SBO_RNA:
                prod_rate = (f"{name_to_symbol['Transcription rate'][0][0:-3]}"
                             f"{{{species}}}")
            elif f_type == sbol3.SBO_PROTEIN:
                prod_rate = (f"{name_to_symbol['Transcription-translation rate'][0][0:-3]}"
                             f"{{{species}}}")
            else:
                raise ValueError(f'Cannot handle type {tyto.SBO.get_term_by_uri(f_type)} in {feature_participation[0]}')
            return f'+ {prod_rate}{modulation}{context}'
        else:
            logging.warning(f'Cannot serialize role in {interaction.identity} of type {tyto.SBO.get_term_by_uri(i_type)}')
    elif i_type == tyto.SBO.cleavage:
        if interaction.name == 'Cas cleavage':
            reactants = [maybe_concentration(f) for f in all_in_role(interaction, sbol3.SBO_REACTANT)]
            rate = '\\casCutRate{{}}'
            if role == sbol3.SBO_REACTANT:
                sign = '-'
            elif role == sbol3.SBO_PRODUCT:
                sign = '+'
            else:
                raise ValueError(f'Unexpected role in {interaction.identity}: {tyto.SBO.get_term_by_uri(role)}')
            return f'{sign} {rate}' + ''.join(reactants)
        else:
            raise ValueError(f'No model for cleavage {interaction.name} in {interaction.identity}')
    elif i_type == sbol3.SBO_DEGRADATION:
        if len(interaction.participations) != 1:
            raise ValueError(f'Degradation assumed to have 1 participant, '
                             f'found {len(interaction.participations)} in '
                             f'{interaction.identity}')
        if f_type == sbol3.SBO_RNA:
            deg_rate = name_to_symbol['gRNA degradation'][0]
        elif f_type == sbol3.SBO_PROTEIN or f_type == sbol3.SBO_DNA or f_type == sbol3.SBO_NON_COVALENT_COMPLEX:
            deg_rate = name_to_symbol['Stable Molecule Dilution'][0]
        else:
            deg_rate = name_to_symbol[interaction.name][0]
        return f'- {deg_rate}{maybe_concentration(feature)}'
    elif i_type == sbol3.SBO_NON_COVALENT_BINDING:
        reactants = [maybe_concentration(f) for f in all_in_role(interaction, sbol3.SBO_REACTANT)]
        rate = name_to_symbol[interaction.name][0] # TODO: move this into actual parameters rather than name
        if role == sbol3.SBO_REACTANT:
            sign = '-'
        elif role == sbol3.SBO_PRODUCT:
            sign = '+'
        else:
            raise ValueError(f'Cannot handle type {tyto.SBO.get_term_by_uri(f_type)} in {interaction.identity}')
        return f'{sign} {rate}' + ''.join(reactants)
    elif i_type == sbol3.SBO_INHIBITION or i_type == sbol3.SBO_STIMULATION:
        # Pass for the regulation interactions that are taken care of in the other function, so you don't get a warning
        pass
    elif i_type == RECOMBINATION:
        if role == sbol3.SBO_MODIFIER or role == sbol3.SBO_MODIFIED:
            return None  # no effect on Cre concentration, not modeling excised element
        reactant = in_role(interaction, sbol3.SBO_REACTANT)
        recombinase = maybe_concentration(in_role(interaction, sbol3.SBO_MODIFIER))
        ct = containers[reactant]
        if len(ct) != 1:
            raise ValueError(f'Recombination expected 1 context, got {len(ct)} in {interaction.identity}')
        context = ct[0]
        rate = name_to_symbol[interaction.name][0] # TODO: move this into actual parameters rather than name
        if role == sbol3.SBO_REACTANT:
            sign = '-'
        elif role == sbol3.SBO_PRODUCT:
            sign = '+'
        else:
            raise ValueError(f'Cannot handle type {tyto.SBO.get_term_by_uri(f_type)} in {interaction.identity}')
        return f'{sign} {rate} {maybe_concentration(reactant)} {recombinase}^4 + ' \
               f'\\frac{{{maybe_concentration(feature)}}}{{{maybe_concentration(context)}}} {differential(context)}'
    else:
        logging.warning(f'Cannot serialize interaction {interaction.identity} of type {tyto.SBO.get_term_by_uri(i_type)}')
        return None


def make_latex_model(system: sbol3.Component) -> str:
    """Generate a set of LaTeX equations for the identified system:

    :param system: system for which a model is to be generated
    :return: string serialization of LaTeX equation collection
    """
    # for each feature, collect all of the interactions and constraints that it
    # participates in
    interactions = {f: [i for i in system.interactions
                        if [p for p in i.participations
                            if p.participant == f.identity]]
                    for f in system.features}
    direct_containers = {f: [c.subject.lookup() for c in system.constraints
                             if c.restriction == sbol3.SBOL_CONTAINS
                             and c.object == f.identity]
                         for f in system.features}
    containers = transitive_closure(direct_containers)
    regulators = {f: [c.subject.lookup() for c in system.constraints
                      if c.restriction == sbol3.SBOL_MEETS
                      and c.object == f.identity]
                  for f in system.features}
    regulation = {f: list(itertools.chain(*(interactions[r] for r in regulators[f]))) for f in regulators}

    # generate an ODE based on the roles in the interactions
    equation_latex = []
    for f in id_sort(system.features):
        interaction_terms = [t for t in [interaction_to_term(f, i, interactions, regulation, containers) for i in id_sort(interactions[f])] if t]
        # If there is at least one term, then add an equation
        if interaction_terms:
            equation_latex.append(f'{differential(f)} & = ' + ' '.join(sorted(interaction_terms)).removeprefix('+'))

    ## Generate the actual document
    # write section header
    latex =  f'\\subsection{{{system.name or system.display_id}}}\n\\label{{s:{system.display_id}}}\n'
    latex += f'% Equations generated from {system.identity}\n\n'
    if system.description:
        latex += f'{system.description}\n\n'
    # write equations
    latex += f'\\begin{{align}}\n'
    latex += '\\\\\n'.join(equation_latex)
    latex += f'\n\\end{{align}}\n\n'

    return latex
