import numpy as np
import re

def print_name(model_name):
    """Print the model name."""
    print('"name": "{}",'.format(model_name))

def take_number_from_string(main_string, substring):
    """Count the elements from a string."""
    if substring == main_string:
        return 1
    elif substring in main_string:
        return int(main_string.replace(substring, ''))
    else:
        raise ValueError("Check the chemical formula.")

def count_sites(formula):
    """Count the number of sites in the chemical formula."""
    open_parentheses = formula.count('(')
    close_parentheses = formula.count(')')
    if open_parentheses == close_parentheses:
        return open_parentheses
    else:
        raise ValueError("Check the chemical formula.")

def count_elems_in_formula(formula, nsites, elems, comp):
    """Find substrings between parentheses and count the elemets in the formula."""
    substrings = re.findall('\((.*?)\)', formula)
    for i in range(nsites):
        substring = substrings[i].split()
        for j in range(len(substring)): 
            for k in range(len(elems)):
                if elems[k] in substring[j]:
                    nelem = take_number_from_string(substring[j], elems[k])
                    comp[i, k] += nelem
    return comp

def print_formula(key, comp, flag):
    """Print the chemical formula."""
    print('    "{}": ['.format(key), end="")
    for i in range(np.shape(comp)[0]):
        print("[", end="")
        for j in range(np.shape(comp)[1]):
            if j < np.shape(comp)[1]-1:
                print(comp[i, j], end=", ")
            else:
                if i < np.shape(comp)[0]-1:
                    print(comp[i, j], end="], ")
                else:
                    print(comp[i, j], end="]")
    if flag:
        print("]", end=",\n")
    else:
        print("]", end="\n")
    

def print_comps(formulas, elems, flags):
    
    flag = True
    aux = 0
    for key, formula in formulas.items():
        
        # Make all the strings upppercase to avoid errors
        for e in range(len(elems)):
            elems[e] = elems[e].upper()
        formula = formula.upper()

        # Count the number of sites
        nsites = count_sites(formula)

        if aux == 0:
            print('"sites": {},'.format(nsites))
            print('"endmembers": {')
        if aux == len(formulas)-1:
            flag = False

        # Initialize the composition matrix
        comp = np.zeros((nsites, len(elems)))

        comp = count_elems_in_formula(formula, nsites, elems, comp)
        comp = comp*flags
        print_formula(key, comp, flag)
        aux += 1
    print('},')

def print_margules(margules, endmemers):
    """Print the margules parameters as vectors."""
    print('"margules": {')
    aux = 1
    for i in range(np.shape(margules)[0]):
        for j in range(aux, np.shape(margules)[1]):
            if i < np.shape(margules)[0]-2:
                str = '    "{},{}": {:.3e},'.format(endmemers[i], endmemers[j], margules[i][j])
            else:
                str = '    "{},{}": {:.3e}'.format(endmemers[i], endmemers[j], margules[i][j])
            print(str)
        aux += 1
    print('}')

def change_margules(margules, endmemers, endmember_a, endmember_b, value):
    """Change the margules parameters."""
    i = endmemers.index(endmember_a)
    j = endmemers.index(endmember_b)
    margules[i][j] = value
    margules[j][i] = value

def print_margules_burnman(margules, endmemers):
    """Print the margules parameters as vectors."""
    print('margules_burnman: [', end="")
    aux = 1
    for i in range(np.shape(margules)[0]-1):
        print("[", end="")
        for j in range(aux, np.shape(margules)[1]):
            
            if i < np.shape(margules)[0]-2 and j < np.shape(margules)[1]-1:
                str = '{:.3e}, '.format(margules[i][j])
            else:
                str = '{:.3e}'.format(margules[i][j])
            print(str, end="")
        aux += 1
        if i < np.shape(margules)[0]-2:
            print("], ", end="")
        else:
            print("]", end="")
    print(']')