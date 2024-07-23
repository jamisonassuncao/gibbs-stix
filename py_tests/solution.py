from formula import *

# List of elements to count
elems = ["Si", "Ca", "Al", "Fe", "Mg", "Na"]
# Flags to multiply the <elems> (effectively removing some elements that should not be counted)
# flags = [0.00, 1.00, 0.00, 1.00, 1.00, 1.00] 
flags = [1.00, 1.00, 1.00, 1.00, 1.00, 1.00] 

# ################################################################################

# # Model name
# model_name = "plagioclase"

# # Chemical formulas dictionary
# formulas = {
#     "anorthite":"(Ca)(Al2)Si2O8",
#     "albite":"(Na)(Al Si3)O8",
# }

# # Specify the endmembers
# endmemers = list(formulas.keys())

# # Initialize the margules parameters
# margules = np.zeros((len(endmemers), len(endmemers)))

# # Change the margules parameters individualy
# change_margules(margules, endmemers, "anorthite", "albite", 26.0e3)

################################################################################



# # Model name
# model_name = "orthopyroxene"

# # Chemical formulas dictionary
# formulas = {
#     "enstatite":"(Mg)(Mg)Si2O6",
#     "ferrosilite":"(Fe)(Fe)Si2O6",
#     "mg-tschermak":"(Mg)(Al)SiAlO6",
#     "ortho-diopside":"(Ca)(Mg)Si2O6",
# }

# # Specify the endmembers
# endmemers = list(formulas.keys())

# # Initialize the margules parameters
# margules = np.zeros((len(endmemers), len(endmemers)))

# # Change the margules parameters individualy
# change_margules(margules, endmemers, "ortho-diopside", "mg-tschermak", 48.0e3)
# change_margules(margules, endmemers, "ortho-diopside", "enstatite", 32.1e3)

################################################################################

# # Model name
# model_name = "clinopyroxene"

# # Chemical formulas dictionary
# formulas = {
#     "diopside": "(Ca)(Mg)(Si2)O6",
#     "hedenbergite": "(Ca)(Fe)(Si2)O6",
#     "clinoenstatite": "(Mg)(Mg)(Si2)O6",
#     "ca-tschermak": "(Ca)(Al)(Si Al)O6",
#     "jadeite": "(Na)(Al)(Si2)O6",
# }

# # Specify the endmembers
# endmemers = list(formulas.keys())

# # Initialize the margules parameters
# margules = np.zeros((len(endmemers), len(endmemers)))

# # Change the margules parameters individualy
# change_margules(margules, endmemers, "clinoenstatite", "diopside", 24.7e3)
# change_margules(margules, endmemers, "jadeite", "diopside", 24.3e3)
# change_margules(margules, endmemers, "ca-tschermak", "diopside", 26e3)
# change_margules(margules, endmemers, "clinoenstatite", "ca-tschermak", 60.6e3)
# change_margules(margules, endmemers, "clinoenstatite", "hedenbergite", 24.7e3)
# change_margules(margules, endmemers, "jadeite", "ca-tschermak", 10.0e3)

# ################################################################################

# # Model name
# model_name = "hp-clinopyroxene"

# # Chemical formulas dictionary
# formulas = {
#     "hp-clinoenstatite": "(Mg2)Si2O6",
#     "hp-clinoferrosilite": "(Fe2)Si2O6",
# }

# # Specify the endmembers
# endmemers = list(formulas.keys())

# # Initialize the margules parameters
# margules = np.zeros((len(endmemers), len(endmemers)))

# # Change the margules parameters individualy
# change_margules(margules, endmemers, "hp-clinoenstatite", "hp-clinoferrosilite", 0.0e3)

# ################################################################################

# # Model name
# model_name = "akimotoite"

# # Chemical formulas dictionary
# formulas = {
#     "mg-akimotoite": "(Mg)(Si)O3",
#     "fe-akimotoite": "(Fe)(Si)O3",
#     "corundum": "(Al)(Al)O3",
# }

# # Specify the endmembers
# endmemers = list(formulas.keys())

# # Initialize the margules parameters
# margules = np.zeros((len(endmemers), len(endmemers)))

# # Change the margules parameters individualy
# change_margules(margules, endmemers, "mg-akimotoite", "corundum", 66.0e3)

################################################################################

# Model name
model_name = "garnet"

# Chemical formulas dictionary
formulas = {
    "pyrope": "(Mg3)(Al)(Al)Si3O12",
    "almandine": "(Fe3)(Al)(Al)Si3O12",
    "grossular": "(Ca3)(Al)(Al)Si3O12",
    "mg-majorite": "(Mg3)(Mg)(Si)Si3O12",
    "jd-majorite": "(Na2 Al)(Al)(Si)Si3O12",
}

# Specify the endmembers
endmemers = list(formulas.keys())

# Initialize the margules parameters
margules = np.zeros((len(endmemers), len(endmemers)))

# Change the margules parameters individualy
change_margules(margules, endmemers, "grossular", "mg-majorite", 58e3)
change_margules(margules, endmemers, "grossular", "pyrope", 30e3)
change_margules(margules, endmemers, "pyrope", "mg-majorite", 21.3e3)

################################################################################

# # Model name
# model_name = "perovskite"

# # Chemical formulas dictionary
# formulas = {
#     "mg-perovskite": "(Mg)(Si)O3",
#     "fe-perovskite": "(Fe)(Si)O3",
#     "al-perovskite": "(Al)(Al)O3",
# }

# # Specify the endmembers
# endmemers = list(formulas.keys())

# # Initialize the margules parameters
# margules = np.zeros((len(endmemers), len(endmemers)))

# # Change the margules parameters individualy
# change_margules(margules, endmemers, "mg-perovskite", "al-perovskite", 116e3)

# ################################################################################

# # Model name
# model_name = "post-perovskite"

# # Chemical formulas dictionary
# formulas = {
#     "mg-post-perovskite": "(Mg)(Si)O3",
#     "fe-post-perovskite": "(Fe)(Si)O3",
#     "al-post-perovskite": "(Al)(Al)O3",
# }

# # Specify the endmembers
# endmemers = list(formulas.keys())

# # Initialize the margules parameters
# margules = np.zeros((len(endmemers), len(endmemers)))

# # Change the margules parameters individualy
# change_margules(margules, endmemers, "mg-post-perovskite", "al-post-perovskite", 60.0e3)

# ################################################################################

# # Model name
# model_name = "magneto-wustite"

# # Chemical formulas dictionary
# formulas = {
#     "periclase": "(Mg)O",
#     "wustite": "(Fe)O",
# }

# # Specify the endmembers
# endmemers = list(formulas.keys())

# # Initialize the margules parameters
# margules = np.zeros((len(endmemers), len(endmemers)))

# # Change the margules parameters individualy
# change_margules(margules, endmemers, "periclase", "wustite", 13e3)
################################################################################

# # Model name
# model_name = "magneto-wustite"

# # Chemical formulas dictionary
# formulas = {
#     "mg-ca-ferrite": "(Mg)(Al)AlO4",
#     "fe-ca-ferrite": "(Fe)(Al)AlO4",
#     "na-ca-ferrite": "(Na)(Si)AlO4",
# }

# # Specify the endmembers
# endmemers = list(formulas.keys())

# # Initialize the margules parameters
# margules = np.zeros((len(endmemers), len(endmemers)))

# # Change the margules parameters individualy
# change_margules(margules, endmemers, "mg-ca-ferrite", "fe-ca-ferrite", 0.0)

################################################################################
print_name(model_name)
print_comps(formulas, elems, flags)
print_margules(margules, endmemers)
print_margules_burnman(margules, endmemers)