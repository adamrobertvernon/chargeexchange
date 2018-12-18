import numpy as np
import pandas as pd
import re #get rid of nun numeric strings

def fill_zeros_with_last(arr):
	last_val = None # I don't really care about the initial value
	for i in range(arr.size):
		if arr[i] and arr[i]!='':
			last_val = arr[i]
		else:
			arr[i] = last_val

def load_levels_lines(B_string, second_ce, dist, database, ele_B, charge_state, I_A):

  if second_ce:
    B_string=ele_sym_B+'I'

  I_levels_file="datatables/"+B_string+'cm-1.nist'

  I_levels = open(I_levels_file,'rb')

  if dist!=0:
    if database == 'nist':
      I_lines_file="datatables/"+B_string+'linescm-1.nist'
    if database == 'kurucz':
      I_lines_file="datatables/"+B_string+'Ilinescm-1.kurucz'
    try:
      I_lines = open(I_lines_file,'rb')
    except Exception as e:
      print(e)


  I_B=ele_B.ionenergies[charge_state]

  if I_A == I_B: #resonant charge exchange not treated
      I_B=I_B+0.000005

  print("I_A:",I_A,"I_B:",I_B)


  data = np.genfromtxt(I_levels,delimiter = '|',dtype=str,skip_header=4,autostrip=1)
  data=data[np.array([len(entry) > 1 for entry in data.T[3]])]# get rid of length 0 strings i.e. empty energies

  atom_config = list(filter(None,data.T[0]))
  term = np.array(list(data.T[1]))#n.b. will have whitespace
  # term = list(data.T[1])#n.b. will have whitespace
  J = [entry.replace("[","").replace("]","").replace(" ","").replace("?","").replace("x","").replace("+","") for entry in list(data.T[2])]
  level_raw = [entry.replace("[","").replace("]","").replace(" ","").replace("?","").replace("x","").replace("+","") for entry in list(data.T[3])]
  fill_zeros_with_last(term)


  if dist!=0:
    if "nist" in I_lines_file:
      print("loading transitions...")
      data = np.genfromtxt(I_lines,delimiter = '|',dtype=str,skip_header=6,autostrip=1)
      print("done loading")
      # print(data)
      A_ki = list(filter(None,data.T[3])) #s^-1
      # lines = list(filter(None,data.T[5])) #s^-1
      lines = [entry.replace("[","").replace("]","").replace(" ","").replace("?","").replace("x","").replace("+","") for entry in list(filter(None,data.T[5]))] # s-1

      E_i=[]
      E_k=[]

      for x in lines:
        transition=x.split("-")
        E_i.append(float(transition[0]))
        E_k_val=transition[1].replace("[","").replace("]","").replace(" ","")
        # print(E_k_val)
        # E_k_val=re.sub("[^0-9]", "", transition[1])
        E_k.append(float(E_k_val))

    elif "kurucz" in I_lines_file:
      data = np.genfromtxt(I_lines,dtype=str,skip_header=4,autostrip=1) #default delimiter is strings of whitespace
      A_ki = list(filter(None,data.T[2])) #s^-1
      E_i=[float(i) for i in list(filter(None,data.T[6]))]
      E_k=[float(i) for i in list(filter(None,data.T[8]))]
      # lines=list(zip(E_i,E_k))
    else:
      print("not .nist or .kurucz??!")


  ####fills in the blanks terms and creates Ls and Ss####
  def spec_to_L(spec):

    spec_dict={"S":0,"P":1,"D":2,"F":3,"G":4,"H":5,"I":6,"K":7,"L":8,"M":9}

    for letter in spec:
      if letter in spec_dict.keys():
        return spec_dict[letter]
    else:
      print("Spec notation not recognised")


    # if spec[0] =='[':
    if "[" in spec:
      val=spec[spec.index("[") + 1:spec.rindex("]")]
      return float(eval(val))
    else:
      return spec

  L, S, level, skipped_level=[], [], [], []

  for index, x in enumerate(term):
    print("term", x)
    # temp_term=x

    if x == '':
      # if J[index] != '':
      # 	term[index]=temp_term
      temp_term="thiswontfloat"
    else:
      temp_term=x

    if ' ' in temp_term: temp_term=temp_term[temp_term.index(' ')+1:] # fixes case e.g. "a 7S*", "c 5G" ...

    if J[index] != '':
      try:
        L_mom=float(spec_to_L(temp_term[1:]))
        L.append(L_mom)
        # print("level_raw[index]", level_raw[index], len(level_raw[index]))
        level.append(level_raw[index])
      except Exception as e: #if ang. mom. assignments unknown
        print(e)
        print("skipping level", temp_term)
        if level_raw[index][0] =='[':
          val=level_raw[index][level_raw[index].index("[") + 1:level_raw[index].rindex("]")]
          skipped_level.append(float(val))
        else:
          # val=re.sub("[^0-9]", "", level_raw[index])
          # skipped_level.append(float(val))
          skipped_level.append(float(level_raw[index]))
      try:
        S.append((int(temp_term[0])-1)/2)
      except:
        pass
    elif level_raw[index] != '':
      term[index-1]=temp_term
      try:
        L_mom=float(spec_to_L(temp_term[1:]))
        L.append(L_mom)
        level.append(level_raw[index])
      except Exception as e: #if ang. mom. assignments unknown
        print(e)
        print("skipping level", temp_term)
      try:
        S.append((int(temp_term[0])-1)/2)
      except:
        pass


  #removes leftover blanks
  term=list(filter(None,term))
  J=list(filter(None,J))
  ####fills in the blanks terms and creates Ls and Ss####

  if dist!=0:
    E_k_probs=[]#decay out of
    for index_x,x in enumerate(E_k,start=0):
      for index_y,y in enumerate(E_k,start=0):
        try:
          if float(y)==float(x) and not (y in E_k_lines[0]):
            E_k_probs.append([x,E_i[index_y],float(A_ki[index_y])])
        except: #excepts on first time
          if float(y)==float(x):
            E_k_probs.append([x,E_i[index_y],float(A_ki[index_y])])
            #sum_A_ki.append[x,A_ki[index_y]]

  if dist!=0:
    return  I_B, term, J, L, S, level, skipped_level,  E_k_probs
  else:
    return  I_B, term, J, L, S, level, skipped_level
