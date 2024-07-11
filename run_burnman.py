from burnman import Solution, minerals
from burnman.classes.solutionmodel import SymmetricRegularSolution, AsymmetricRegularSolution
from burnman.tools.chemistry import formula_to_string
import numpy as np

def fml(mineral):
    str = formula_to_string(mineral.formula)
    return str

TO_PA = 100_000.0
pressure = 1000.0 * TO_PA # bar 2 Pa
temperature = 1000.0 # K

# mineral_a = minerals.SLB_2011.albite()
# mineral_b = minerals.SLB_2011.anorthite()
# model = Solution(name = 'Plagioclase', solution_model = SymmetricRegularSolution(endmembers= [[mineral_a, fml(mineral_a)], [mineral_b, fml(mineral_b)]], energy_interaction=[[26.0e3]]))
# model.set_composition([0.5, 0.5])
# model.set_state(pressure, temperature)
# mineral_a.set_state(pressure, temperature)
# mineral_b.set_state(pressure, temperature)

<<<<<<< HEAD
mineral_a = minerals.SLB_2011.spinel()
mineral_b = minerals.SLB_2011.hercynite()
model = Solution(name = 'Spinel', solution_model = SymmetricRegularSolution(endmembers = [[mineral_a, '[Mg3Al1][Al7Mg1]1O16'], [mineral_b, '[Fe3Al1]1[Al7Fe1]1O16']], energy_interaction=[[5000.0]]))
model.set_composition([0.5, 0.5])
model.set_state(pressure, temperature)
mineral_a.set_state(pressure, temperature)
mineral_b.set_state(pressure, temperature)
=======
# mineral_a = minerals.SLB_2011.spinel()
# mineral_b = minerals.SLB_2011.hercynite()
# model = Solution(name = 'Spinel', solution_model = SymmetricRegularSolution(endmembers = [[mineral_a, '[Mg3Al1][Al7Mg1]1O16'], [mineral_b, '[Fe3Al1]1[Al7Fe1]1O16']], energy_interaction=[[5000.0]]))
# model.set_composition([0.5, 0.5])
# model.set_state(pressure, temperature)
# mineral_a.set_state(pressure, temperature)
# mineral_b.set_state(pressure, temperature)
>>>>>>> 810b2137dfaab5d6d63e10c27f821e69d2f5324b

# mineral_a = minerals.SLB_2011.forsterite()
# mineral_b = minerals.SLB_2011.fayalite()
# model = Solution(name = 'Olivine', solution_model = SymmetricRegularSolution(endmembers= [[mineral_a,'[Mg2]Si1O4'], [mineral_b, '[Fe2]Si1O4']], energy_interaction=[[7.6e3]]))
<<<<<<< HEAD
# model.set_composition([0.7, 0.3])
=======
# model.set_composition([0.5, 0.5])
>>>>>>> 810b2137dfaab5d6d63e10c27f821e69d2f5324b
# model.set_state(pressure, temperature)
# mineral_a.set_state(pressure, temperature)
# mineral_b.set_state(pressure, temperature)

# mineral_a = minerals.SLB_2011.mg_wadsleyite()
# mineral_b = minerals.SLB_2011.fe_wadsleyite()
# model = Solution(name = 'Wadsleyite', solution_model = SymmetricRegularSolution(endmembers= [[mineral_a,'[Mg2]Si1O4'], [mineral_b, '[Fe2]Si1O4']], energy_interaction=[[16.5e3]]))
# model.set_composition([0.5, 0.5])
# model.set_state(pressure, temperature)
# mineral_a.set_state(pressure, temperature)
# mineral_b.set_state(pressure, temperature)

# mineral_a = minerals.SLB_2011.mg_ringwoodite()
# mineral_b = minerals.SLB_2011.fe_ringwoodite()
# model = Solution(name = 'Ringwoodite', solution_model = SymmetricRegularSolution(endmembers= [[mineral_a,'[Mg2]Si1O4'], [mineral_b, '[Fe2]Si1O4']], energy_interaction=[[9.1e3]]))
# model.set_composition([0.5, 0.5])
# model.set_state(pressure, temperature)
# mineral_a.set_state(pressure, temperature)
# mineral_b.set_state(pressure, temperature)

# mineral_a = minerals.SLB_2011.enstatite()
# mineral_b = minerals.SLB_2011.ferrosilite()
# mineral_c = minerals.SLB_2011.mg_tschermaks()
# mineral_d = minerals.SLB_2011.ortho_diopside()
# model = Solution(name = 'Orthopyroxene', solution_model = SymmetricRegularSolution(endmembers = [[mineral_a,'[Mg_1][Mg_1]Si_2O_6'], [mineral_b, '[Fe_1][Fe_1]Si_2O_6'], [mineral_c, '[Mg_1][Al_1]Si_1Al_1O_6'], [mineral_d, '[Ca_1][Mg_1]Si_2O_6']], energy_interaction=[[0.000e+00, 0.000e+00, 3.210e+04], [0.000e+00, 0.000e+00], [4.800e+04]]))
# model.set_composition([0.25, 0.25, 0.25, 0.25])
# model.set_state(pressure, temperature)
# mineral_a.set_state(pressure, temperature)
# mineral_b.set_state(pressure, temperature)
# mineral_c.set_state(pressure, temperature)
# mineral_d.set_state(pressure, temperature)

<<<<<<< HEAD
# mineral_a = minerals.SLB_2011.diopside()
# mineral_b = minerals.SLB_2011.hedenbergite()
# mineral_c = minerals.SLB_2011.clinoenstatite()
# mineral_d = minerals.SLB_2011.ca_tschermaks()
# mineral_e = minerals.SLB_2011.jadeite()
# model = Solution(name= 'Clinopyroxene', solution_model = AsymmetricRegularSolution(endmembers= [[mineral_a, '(Ca)(Mg)(Si2)O6'], [mineral_b, '(Ca)(Fe)(Si2)O6'], [mineral_c, '(Mg)(Mg)(Si2)O6'], [mineral_d, '(Ca)(Al)(Si Al)O6'], [mineral_e, '(Na)(Al)(Si2)O6']], alphas=[1.0, 1.0, 1.0, 3.5, 1.0], energy_interaction=[[0.000e+00, 2.470e+04, 2.600e+04, 2.430e+04],[2.470e+04, 0.000e+00, 0.000e+00],[6.060e+04, 0.000e+00],[1.000e+04]]))
# model.set_composition([0.2, 0.2, 0.2, 0.2, 0.2])
# model.set_state(pressure, temperature)
# mineral_a.set_state(pressure, temperature)
# mineral_b.set_state(pressure, temperature)
# mineral_c.set_state(pressure, temperature)
# mineral_d.set_state(pressure, temperature)
# mineral_e.set_state(pressure, temperature)
=======
mineral_a = minerals.SLB_2011.diopside()
mineral_b = minerals.SLB_2011.hedenbergite()
mineral_c = minerals.SLB_2011.clinoenstatite()
mineral_d = minerals.SLB_2011.ca_tschermaks()
mineral_e = minerals.SLB_2011.jadeite()
model = Solution(name= 'Clinopyroxene', solution_model = AsymmetricRegularSolution(endmembers= [[mineral_a, '(Ca)(Mg)(Si2)O6'], [mineral_b, '(Ca)(Fe)(Si2)O6'], [mineral_c, '(Mg)(Mg)(Si2)O6'], [mineral_d, '(Ca)(Al)(Si Al)O6'], [mineral_e, '(Na)(Al)(Si2)O6']], alphas=[1.0, 1.0, 1.0, 3.5, 1.0], energy_interaction=[[0.000e+00, 2.470e+04, 2.600e+04, 2.430e+04],[2.470e+04, 0.000e+00, 0.000e+00],[6.060e+04, 0.000e+00],[1.000e+04]]))
model.set_composition([0.2, 0.2, 0.2, 0.2, 0.2])
model.set_state(pressure, temperature)
mineral_a.set_state(pressure, temperature)
mineral_b.set_state(pressure, temperature)
mineral_c.set_state(pressure, temperature)
mineral_d.set_state(pressure, temperature)
mineral_e.set_state(pressure, temperature)
>>>>>>> 810b2137dfaab5d6d63e10c27f821e69d2f5324b

# mineral_a = minerals.SLB_2011.hp_clinoenstatite()
# mineral_b = minerals.SLB_2011.hp_clinoferrosilite()
# model = Solution(name = 'Ringwoodite', solution_model = SymmetricRegularSolution(endmembers= [[mineral_a,'[Mg2]Si2O6'], [mineral_b, '[Fe2]Si2O6']], energy_interaction=[[0.0e3]]))
# model.set_composition([0.5, 0.5])
# model.set_state(pressure, temperature)
# mineral_a.set_state(pressure, temperature)
# mineral_b.set_state(pressure, temperature)

# mineral_a = minerals.SLB_2011.ca_perovskite()
# model = Solution(name = 'Ca-Perovskite', solution_model = SymmetricRegularSolution(endmembers= [[mineral_a,'CaSiO3']], energy_interaction=[[0.0e3]]))
# model.set_composition([1.0])
# model.set_state(pressure, temperature)
# mineral_a.set_state(pressure, temperature)

# mineral_a = minerals.SLB_2011.quartz()
# model = Solution(name = 'Quartz', solution_model = SymmetricRegularSolution(endmembers= [[mineral_a,fml(mineral_a)]], energy_interaction=[[0.0]]))
# model.set_composition([1.0])
# model.set_state(pressure, temperature)
# mineral_a.set_state(pressure, temperature)

# mineral_a = minerals.SLB_2011.coesite()
# model = Solution(name = 'Coesite', solution_model = SymmetricRegularSolution(endmembers= [[mineral_a,'SiO2']], energy_interaction=[[0.0e3]]))
# model.set_composition([1.0])
# model.set_state(pressure, temperature)
# mineral_a.set_state(pressure, temperature)

# mineral_a = minerals.SLB_2011.stishovite()
# model = Solution(name = 'Stishovite', solution_model = SymmetricRegularSolution(endmembers= [[mineral_a,'SiO2']], energy_interaction=[[0.0e3]]))
# model.set_composition([1.0])
# model.set_state(pressure, temperature)
# mineral_a.set_state(pressure, temperature)

# mineral_a = minerals.SLB_2011.seifertite()
# model = Solution(name = 'Seifertite', solution_model = SymmetricRegularSolution(endmembers= [[mineral_a, formula_to_string(mineral_a.formula)]], energy_interaction=[[0.0e3]]))
# model.set_composition([1.0])
# model.set_state(pressure, temperature)
# mineral_a.set_state(pressure, temperature)

# mineral_a = minerals.SLB_2011.mg_perovskite()
# mineral_b = minerals.SLB_2011.fe_perovskite()
# mineral_c = minerals.SLB_2011.al_perovskite()
# model = Solution(name = 'Perovskite', 
#                  solution_model = AsymmetricRegularSolution(endmembers= [[mineral_a, fml(mineral_a)], 
#                                                                         [mineral_b, fml(mineral_b)], 
#                                                                         [mineral_c, fml(mineral_c)]], 
#                                                            alphas = [1.00, 1.00, 0.39],
#                                                            energy_interaction=[[0.000e+00, 1.160e+05], 
#                                                                                [0.000e+00]]))
# model.set_composition([0.3, 0.3, 0.4])
# model.set_state(pressure, temperature)
# mineral_a.set_state(pressure, temperature)
# mineral_b.set_state(pressure, temperature)
# mineral_c.set_state(pressure, temperature)

# mineral_a = minerals.SLB_2011.mg_post_perovskite()
# mineral_b = minerals.SLB_2011.fe_post_perovskite()
# mineral_c = minerals.SLB_2011.al_post_perovskite()
# model = Solution(name = 'Perovskite', 
#                  solution_model = SymmetricRegularSolution(endmembers= [[mineral_a, fml(mineral_a)], 
#                                                                         [mineral_b, fml(mineral_b)], 
#                                                                         [mineral_c, fml(mineral_c)]], 
#                                                            energy_interaction=[[0.000e+00, 6.000e+04], [0.000e+00]]))
# model.set_composition([0.3, 0.3, 0.4])
# model.set_state(pressure, temperature)
# mineral_a.set_state(pressure, temperature)

# mineral_a = minerals.SLB_2011.periclase()
# mineral_b = minerals.SLB_2011.wuestite()
# model = Solution(name = 'Magnesio-wustite', 
#                  solution_model = SymmetricRegularSolution(endmembers= [[mineral_a, fml(mineral_a)], 
#                                                                         [mineral_b, fml(mineral_b)]], 
#                                                            energy_interaction=[[1.300e+04]]))
# model.set_composition([0.5, 0.5])
# model.set_state(pressure, temperature)
# mineral_a.set_state(pressure, temperature)
# mineral_b.set_state(pressure, temperature)

# mineral_a = minerals.SLB_2011.mg_ca_ferrite()
# mineral_b = minerals.SLB_2011.fe_ca_ferrite()
# mineral_c = minerals.SLB_2011.na_ca_ferrite()
# model = Solution(name = 'Ca-Ferrite', 
#                  solution_model = SymmetricRegularSolution(endmembers= [[mineral_a, fml(mineral_a)], 
#                                                                         [mineral_b, fml(mineral_b)],
#                                                                         [mineral_c, fml(mineral_c)]], 
#                                                            energy_interaction=[[00e+04, 00e+04],[00e+04]]))
# model.set_composition([0.3, 0.3, 0.4])
# model.set_state(pressure, temperature)
# mineral_a.set_state(pressure, temperature)
# mineral_b.set_state(pressure, temperature)
# mineral_c.set_state(pressure, temperature)

# mineral_a = minerals.SLB_2011.kyanite()
# model = Solution(name = 'Ca-Ferrite', 
#                  solution_model = SymmetricRegularSolution(endmembers= [[mineral_a, fml(mineral_a)]], 
#                                                            energy_interaction=[[00e+04]]))
# model.set_composition([1.0])
# model.set_state(pressure, temperature)
# mineral_a.set_state(pressure, temperature)

# mineral_a = minerals.SLB_2011.nepheline()
# model = Solution(name = 'Ca-Ferrite', 
#                  solution_model = SymmetricRegularSolution(endmembers= [[mineral_a, fml(mineral_a)]], 
#                                                            energy_interaction=[[00e+04]]))
# model.set_composition([1.0])
# model.set_state(pressure, temperature)
# mineral_a.set_state(pressure, temperature)

# mineral_a = minerals.SLB_2011.mg_akimotoite()
# mineral_b = minerals.SLB_2011.fe_akimotoite()
# mineral_c = minerals.SLB_2011.corundum()
# model = Solution(name = 'Akimotoite', 
#                  solution_model = SymmetricRegularSolution(endmembers= [[mineral_a, fml(mineral_a)], 
#                                                                         [mineral_b, fml(mineral_b)],
#                                                                         [mineral_c, fml(mineral_c)]], 
#                                                            energy_interaction=[[0.000e+00, 6.600e+04], [0.000e+00]]))
# model.set_composition([0.3, 0.3, 0.4])
# model.set_state(pressure, temperature)
# mineral_a.set_state(pressure, temperature)
# mineral_b.set_state(pressure, temperature)
# mineral_c.set_state(pressure, temperature)

# mineral_a = minerals.SLB_2011.anorthite()
# mineral_b = minerals.SLB_2011.albite()
# model = Solution(name = 'Plagioclase', 
#                  solution_model = SymmetricRegularSolution(endmembers= [[mineral_a, fml(mineral_a)], 
#                                                                         [mineral_b, fml(mineral_b)]], 
#                                                            energy_interaction=[[26.0e3]]))
# model.set_composition([0.5, 0.5])
# model.set_state(pressure, temperature)
# mineral_a.set_state(pressure, temperature)
# mineral_b.set_state(pressure, temperature)

# mineral_a = minerals.SLB_2011.pyrope()
# mineral_b = minerals.SLB_2011.almandine()
# mineral_c = minerals.SLB_2011.grossular()
# mineral_d = minerals.SLB_2011.mg_majorite()
# mineral_e = minerals.SLB_2011.jd_majorite()
# model = Solution(name = 'Garnet', 
#                  solution_model = SymmetricRegularSolution(endmembers= [[mineral_a, fml(mineral_a)], 
#                                                                         [mineral_b, fml(mineral_b)],
#                                                                         [mineral_c, fml(mineral_c)],
#                                                                         [mineral_d, fml(mineral_d)],
#                                                                         [mineral_e, fml(mineral_e)]], 
#                                                            energy_interaction=[[0.000e+00, 3.000e+04, 2.130e+04, 0.000e+00], [0.000e+00, 0.000e+00, 0.000e+00], [5.800e+04, 0.000e+00], [0.000e+00]]))
# model.set_composition([0.2, 0.2, 0.2, 0.2, 0.2])
# model.set_state(pressure, temperature)
# mineral_a.set_state(pressure, temperature)
# mineral_b.set_state(pressure, temperature)
# mineral_c.set_state(pressure, temperature)
# mineral_d.set_state(pressure, temperature)
# mineral_e.set_state(pressure, temperature)

<<<<<<< HEAD

# S_conf = model.solution_model._configurational_entropy(model.molar_fractions)


print(f'Gibbs: {model.gibbs:,.4f} [J/mol]')
# print(temperature*S_conf)

# formula = "Fe3Al2Si3O12"
#         formula = dictionarize_formula(formula)
#         self.params = {
#             "name": "Almandine",
#             "formula": formula,
#             "equation_of_state": "slb3",
#             "F_0": -4935516.0,
#             "V_0": 0.00011543,
#             "K_0": 1.738963e11,
#             "Kprime_0": 4.91341,
#             "Debye_0": 741.356,
#             "grueneisen_0": 1.06495,
#             "q_0": 1.42169,
#             "G_0": 96000000000.0,
#             "Gprime_0": 1.40927,
#             "eta_s_0": 2.09292,
#             "n": sum(formula.values()),
#             "molar_mass": formula_mass(formula),


=======
print("Gibbs:", model.gibbs)
>>>>>>> 810b2137dfaab5d6d63e10c27f821e69d2f5324b
