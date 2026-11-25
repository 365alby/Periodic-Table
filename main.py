import periodictable as pt
from sympy import Matrix, lcm
import re
import time
import sys

# run .venv/bin/python main.py in venv terminal to run
def molarmass(formula):
    return pt.formula(formula).mass

def molcalc(formula, mass):
    return mass / pt.formula(formula).mass

def parse_formula(formula):
  tokens = re.findall(r'([A-Z][a-z]?|\(|\)|\d+)', formula)
  stack = [{}]
  i = 0

  while i < len(tokens):
      token = tokens[i]

      if token == "(":
          stack.append({})

      elif token == ")":
          group = stack.pop()
          i += 1
          if i < len(tokens) and tokens[i].isdigit():
              mult = int(tokens[i])
          else:
              mult = 1
              i -= 1
          for elem, count in group.items():
              stack[-1][elem] = stack[-1].get(elem, 0) + count * mult

      elif re.match(r'[A-Z][a-z]?', token):
          elem = token
          i += 1
          if i < len(tokens) and tokens[i].isdigit():
              num = int(tokens[i])
          else:
              num = 1
              i -= 1
          stack[-1][elem] = stack[-1].get(elem, 0) + num
      elif token.isdigit():
          pass

      i += 1

  return stack[0]

def balancing(reactants, products):
  compounds = reactants + products
  elements = sorted({elem for c in compounds for elem in parse_formula(c)})
  matrix = []
  for elem in elements:
      row = []
      for c in reactants:
          row.append(parse_formula(c).get(elem, 0))
      for c in products:
          row.append(-parse_formula(c).get(elem, 0))
      matrix.append(row)

  M = Matrix(matrix)
  vec = M.nullspace()[0]
  denominators = [term.as_numer_denom()[1] for term in vec]
  denom_lcm = lcm(denominators)
  coeffs = [int(term * denom_lcm) for term in vec]
  reactant_coeff = coeffs[:len(reactants)]
  product_coeff  = coeffs[len(reactants):]

  return reactant_coeff, product_coeff

def print_balanced(reactants, products, rc, pc):
  def fmt(c, n):
    if n == 1:
        return c
    return f"{n} {c}"
  left = " + ".join(fmt(c, n) for c, n in zip(reactants, rc))
  right = " + ".join(fmt(c, n) for c, n in zip(products, pc))
  return f"{left} -> {right}"
  
def empirical_formula():
    element_count = int(input("Enter the number of elements in the compound: "))

    elements = []
    moles = []

    for i in range(element_count):
        symbol = input(f"Enter element {i+1} symbol: ")
        mass = float(input(f"Enter mass of {symbol} (g): "))

        mol = mass / getattr(pt, symbol).mass
        elements.append(symbol)
        moles.append(mol)

    # divide by smallest
    smallest = min(moles)
    ratios = [m / smallest for m in moles]

    # try multiplying ratios until close to integers
    multiplier = 1
    while multiplier <= 10:
        test = [r * multiplier for r in ratios]
        if all(abs(round(x) - x) < 0.05 for x in test):
            ratios = [round(x) for x in test]
            break
        multiplier += 1

    # build formula
    formula = ""
    for elem, count in zip(elements, ratios):
        if count == 1:
            formula += elem
        else:
            formula += f"{elem}{count}"

    print("Empirical formula:", formula)

while True:
  print("--------------------------------------------")
  print("Welcome to the periodic table calculator!")
  time.sleep(1)
  print("Select your options")
  time.sleep(0.15)
  print("1. Calculate molar mass")
  time.sleep(0.15)
  print("2. Calculate moles")
  time.sleep(0.15)
  print("3. Balance equation")
  time.sleep(0.15)
  print("4. Empirical Formula")
  time.sleep(0.15)
  print("5. Exit")
  time.sleep(0.15)
  choice = input("Enter your choice: ")
  print("--------------------------------------------")
  match choice:
    case "1":
      compound = input("Enter the chemical formula of the compound: ")
      mass = molarmass(compound)
      print(round(mass, 3))
    case "2":
      compound = input("Enter the chemical formula of the compound: ")
      mass = float(input("Enter the mass of the compound: "))
      moles = molcalc(compound, mass)
      print("There are", round(moles, 3), "moles of", compound, "in", mass, "grams.")
    case "3":
      reactants = []
      products = []
      reactantnum = int(input("Enter the number of reactants: "))
      productnum = int(input("Enter the number of products: "))
      for i in range(reactantnum):
        reactant = input(f"Enter the chemical formula of the reactant {i+1}: ")
        reactants.append(reactant)
      for i in range(productnum):
        product = input(f"Enter the chemical formula of the product {i+1}: ")
        products.append(product)
      rc, pc = balancing(reactants, products)
      print(print_balanced(reactants, products, rc, pc))
    case "4":
      empirical_formula()
    case "5":
      for _ in range(1):  # number of cycles
        for dots in ["", ".", "..", "..."]:
            sys.stdout.write("\rExiting the program" + dots + "   ")
            sys.stdout.flush()
            time.sleep(0.5)
      print("\nThank you for using the periodic table calculator!")
      break
    case _:
       print("Invalid choice. Please try again.")