import numpy as np

with open('sto3g.out') as inFile:
	lines = inFile.readlines()

index = 0
while lines[index].find("Basis Dimension") == -1:
	index = index +1

dims = int(lines[index].split()[4])

index = 0
while lines[index].find("groups of distinct atoms") == -1:
	index = index + 1

num_groups = int(lines[index].split()[2])

bs1_dim = 0
bs2_dim = 0
str1 = ""
str2 = ""
index = 0
while lines[index].find(" Group   1") == -1:
	index = index + 1
if num_groups == 1:
	bs1_dim = int(dims/2)
	bs2_dim = int(dims/2)
	str1 = lines[index].split()[8]
else:
	str1 = lines[index].split()[8]
	str2 = lines[index+1].split()[8]

ns1 = ns2 = np1 = np2 = nd1 = nd2 = nf1 = nf2 = 0
ng1 = ng2 = nh1 = nh2 = ni1 = ni2 = 0

sloc = str1.find("s")
if sloc != -1:
	ns1 = int(str1[:sloc])
ploc = str1.find("p")
if ploc != -1:
	np1 = int(str1[sloc+1:ploc])
dloc = str1.find("d")
if dloc != -1:
	nd1 = int(str1[ploc+1:dloc])
floc = str1.find("f")
if floc != -1:
	nf1 = int(str1[dloc+1:floc])
gloc = str1.find("g")
if gloc != -1:
	ng1 = int(str1[floc+1:gloc])
hloc = str1.find("h")
if hloc != -1:
	nh1 = int(str1[gloc+1:hloc])
iloc = str1.find("i")
if iloc != -1:
	ni1 = int(str1[hloc+1:iloc])

if str2 != "":
	sloc2 = str2.find("s")
	if sloc2 != -1:
		ns2 = int(str2[:sloc2])
	ploc2 = str2.find("p")
	if ploc2 != -1:
		np2 = int(str2[sloc2+1:ploc2])
	dloc2 = str2.find("d")
	if dloc2 != -1:
		nd2 = int(str2[ploc2+1:dloc2])
	floc2 = str2.find("f")
	if floc2 != -1:
		nf2 = int(str2[dloc2+1:floc2])
	gloc2 = str2.find("g")
	if gloc2 != -1:
		ng2 = int(str2[floc2+1:gloc2])
	hloc2 = str2.find("h")
	if hloc2 != -1:
		nh2 = int(str2[gloc2+1:hloc2])
	iloc2 = str2.find("i")
	if iloc2 != -1:
		ni2 = int(str2[hloc2+1:iloc2])
	bs1_dim = ns1 + 3*np1 + 5*nd1 + 7*nf1 + 9*ng1 + 11*nh1 + 13*ni1
	bs2_dim = ns2 + 3*np2 + 5*nd2 + 7*nf2 + 9*ng2 + 11*nh2 + 13*ni2
else:
	ns2 = ns1
	np2 = np1
	nd2 = nd1
	nf2 = nf1
	ng2 = ng1
	nh2 = nh1
	ni2 = ni1

ovlp = np.zeros((bs2_dim,bs1_dim))

index = 0
while lines[index].find("OVERLAP MATRIX") == -1:
	index = index + 1
index = index + 2

print("     ", end="")
for i in range(bs1_dim):
	if i < ns1:
		print("      S", end="")
	elif i < ns1+3*np1:
		print("      P", end="")
	elif i < ns1+3*np1+5*nd1:
		print("      D", end="")
	elif i < ns1+3*np1+5*nd1+7*nf1:
		print("      F", end="")
print("")

if (bs1_dim%6 != 0):
	for j in range(bs1_dim//6+1):
		for i in range(bs2_dim):
			if j < bs1_dim // 6:
				ovlp[i,j*6]=float(lines[1+bs1_dim+index+j*(bs1_dim+bs2_dim)+i].split()[1])
				ovlp[i,j*6+1]=float(lines[1+bs1_dim+index+j*(bs1_dim+bs2_dim)+i].split()[2])
				ovlp[i,j*6+2]=float(lines[1+bs1_dim+index+j*(bs1_dim+bs2_dim)+i].split()[3])
				ovlp[i,j*6+3]=float(lines[1+bs1_dim+index+j*(bs1_dim+bs2_dim)+i].split()[4])
				ovlp[i,j*6+4]=float(lines[1+bs1_dim+index+j*(bs1_dim+bs2_dim)+i].split()[5])
				ovlp[i,j*6+5]=float(lines[1+bs1_dim+index+j*(bs1_dim+bs2_dim)+i].split()[6])
			else:
				for k in range(1,bs1_dim%6+1):
					ovlp[i,j*6+k-1]=float(lines[1+bs1_dim+index+j*(bs1_dim+bs2_dim)+i].split()[k])
else:
	for j in range(bs1_dim//6):
		for i in range(bs2_dim):
			if j < bs1_dim // 6:
				ovlp[i,j*6]=float(lines[1+bs1_dim+index+j*(bs1_dim+bs2_dim)+i].split()[1])
				ovlp[i,j*6+1]=float(lines[1+bs1_dim+index+j*(bs1_dim+bs2_dim)+i].split()[2])
				ovlp[i,j*6+2]=float(lines[1+bs1_dim+index+j*(bs1_dim+bs2_dim)+i].split()[3])
				ovlp[i,j*6+3]=float(lines[1+bs1_dim+index+j*(bs1_dim+bs2_dim)+i].split()[4])
				ovlp[i,j*6+4]=float(lines[1+bs1_dim+index+j*(bs1_dim+bs2_dim)+i].split()[5])
				ovlp[i,j*6+5]=float(lines[1+bs1_dim+index+j*(bs1_dim+bs2_dim)+i].split()[6])
			else:
				for k in range(1,bs1_dim%6+1):
					ovlp[i,j*6+k-1]=float(lines[1+bs1_dim+index+j*(bs1_dim+bs2_dim)+i].split()[k])

pindex1 = 0
pindex2 = 0
for i in range(bs2_dim):
	if i < ns2:
		print("    S",end="")
	elif i < ns2+3*np2:
		print("    P",end="")
	elif i < ns2+3*np2+5*nd2:
		print("    D",end="")
	for j in range(bs1_dim):
		if i < ns2 and j < ns1:
			print("{:7.3f}".format(ovlp[i,j]),end="")
		elif i < ns2 and j < ns1+3*np1:
			if pindex1 == 0:
				print("{:7.3f}".format(ovlp[i,j+2]),end="")
				pindex1 = 1
			elif pindex1 == 1:
				print("{:7.3f}".format(ovlp[i,j-1]),end="")
				pindex1 = 2
			else:
				print("{:7.3f}".format(ovlp[i,j-1]),end="")
				pindex1 = 0
		elif i < ns2 + 3*np2 and j < ns1:
			if pindex2 == 0:
				print("{:7.3f}".format(ovlp[i+2,j]),end="")
				if pindex1 == 2:
					pindex2 = 1
			elif pindex2 == 1:
				print("{:7.3f}".format(ovlp[i-1,j]),end="")
				if pindex1 == 2:
					pindex2 = 2
			elif pindex2 == 2:
				print("{:7.3f}".format(ovlp[i-1,j]),end="")
				if pindex1 == 2:
					pindex2 = 0
		elif ( i < ( ns2 + 3*np2 )) and ( j < ( ns1 + 3*np1 )):
			if pindex1 == 0 and pindex2 == 0:
				print("{:7.3f}".format(ovlp[i+2,j+2]),end="")
				pindex1 = 1
			elif pindex1 == 0 and pindex2 == 1:
				print("{:7.3f}".format(ovlp[i-1,j+2]),end="")
				pindex1 = 1
			elif pindex1 == 0 and pindex2 == 2:
				print("{:7.3f}".format(ovlp[i-1,j+2]),end="")
				pindex1 = 1
			elif pindex1 == 1 and pindex2 == 0:
				print("{:7.3f}".format(ovlp[i+2,j-1]),end="")
				pindex1 = 2
			elif pindex1 == 1 and pindex2 == 1:
				print("{:7.3f}".format(ovlp[i-1,j-1]),end="")
				pindex1 = 2
			elif pindex1 == 1 and pindex2 == 2:
				print("{:7.3f}".format(ovlp[i-1,j-1]),end="")
				pindex1 = 2
			elif pindex1 == 2 and pindex2 == 0:
				print("{:7.3f}".format(ovlp[i+2,j-1]),end="")
				pindex1 = 0
				pindex2 = 1
			elif pindex1 == 2 and pindex2 == 1:
				print("{:7.3f}".format(ovlp[i-1,j-1]),end="")
				pindex1 = 0
				pindex2 = 2
			elif pindex1 == 2 and pindex2 == 2:
				print("{:7.3f}".format(ovlp[i-1,j-1]),end="")
				pindex1 = 0
				pindex2 = 0
	print("")

