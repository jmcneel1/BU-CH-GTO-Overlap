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
				ovlp[i,j*6]=float(lines[1+bs1_dim+index+j*(bs1_dim+bs2_dim+1)+i].split()[1])
				ovlp[i,j*6+1]=float(lines[1+bs1_dim+index+j*(bs1_dim+bs2_dim+1)+i].split()[2])
				ovlp[i,j*6+2]=float(lines[1+bs1_dim+index+j*(bs1_dim+bs2_dim+1)+i].split()[3])
				ovlp[i,j*6+3]=float(lines[1+bs1_dim+index+j*(bs1_dim+bs2_dim+1)+i].split()[4])
				ovlp[i,j*6+4]=float(lines[1+bs1_dim+index+j*(bs1_dim+bs2_dim+1)+i].split()[5])
				ovlp[i,j*6+5]=float(lines[1+bs1_dim+index+j*(bs1_dim+bs2_dim+1)+i].split()[6])
			else:
				for k in range(1,bs1_dim%6+1):
					ovlp[i,j*6+k-1]=float(lines[1+bs1_dim+index+j*(bs1_dim+bs2_dim+1)+i].split()[k])
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
dindex1 = 0
dindex2 = 0
findex1 = 0
findex2 = 0
for i in range(bs2_dim):
	if i < ns2:
		print("    S",end="")
	elif i < ns2+3*np2:
		print("    P",end="")
	elif i < ns2+3*np2+5*nd2:
		print("    D",end="")
	elif i < ns2+3*np2+5*nd2+7*nf2:
		print("    F",end="")
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
		elif i < ns2 and j < ns1 + 3*np1 + 5*nd1:
			if dindex1 == 0:
				print("{:7.3f}".format(ovlp[i,j+4]),end="")
				dindex1 = 1
			elif dindex1 == 1:
				print("{:7.3f}".format(ovlp[i,j+1]),end="")
				dindex1 = 2
			elif dindex1 == 2:
				print("{:7.3f}".format(ovlp[i,j-2]),end="")
				dindex1 = 3
			elif dindex1 == 3:
				print("{:7.3f}".format(ovlp[i,j-2]),end="")
				dindex1 = 4
			elif dindex1 == 4:
				print("{:7.3f}".format(ovlp[i,j-1]),end="")
				dindex1 = 0
		elif i < ns2 and j < ns1 + 3*np1 + 5*nd1 + 7*nf1:
			if findex1 == 0:
				print("{:7.3f}".format(ovlp[i,j+6]),end="")
				findex1 = 1
			elif findex1 == 1:
				print("{:7.3f}".format(ovlp[i,j+3]),end="")
				findex1 = 2
			elif findex1 == 2:
				print("{:7.3f}".format(ovlp[i,j]),end="")
				findex1 = 3
			elif findex1 == 3:
				print("{:7.3f}".format(ovlp[i,j-3]),end="")
				findex1 = 4
			elif findex1 == 4:
				print("{:7.3f}".format(ovlp[i,j-3]),end="")
				findex1 = 5
			elif findex1 == 5:
				print("{:7.3f}".format(ovlp[i,j-2]),end="")
				findex1 = 6
			elif findex1 == 6:
				print("{:7.3f}".format(ovlp[i,j-1]),end="")
				findex1 = 0
		elif i < ns2 + 3*np2 and j < ns1:
			if pindex2 == 0:
				print("{:7.3f}".format(ovlp[i+2,j]),end="")
			elif pindex2 == 1:
				print("{:7.3f}".format(ovlp[i-1,j]),end="")
			elif pindex2 == 2:
				print("{:7.3f}".format(ovlp[i-1,j]),end="")
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
		elif i < ns2 + 3*np2 and j < ns1 + 3*np1 + 5*nd1:
			if dindex1 == 0 and pindex2 == 0:
				print("{:7.3f}".format(ovlp[i+2,j+4]),end="")
				dindex1 = 1
			elif dindex1 == 1 and pindex2 == 0:
				print("{:7.3f}".format(ovlp[i+2,j+1]),end="")
				dindex1 = 2
			elif dindex1 == 2 and pindex2 == 0:
				print("{:7.3f}".format(ovlp[i+2,j-2]),end="")
				dindex1 = 3
			elif dindex1 == 3 and pindex2 == 0:
				print("{:7.3f}".format(ovlp[i+2,j-2]),end="")
				dindex1 = 4
			elif dindex1 == 4 and pindex2 == 0:
				print("{:7.3f}".format(ovlp[i+2,j-1]),end="")
				dindex1 = 0
			elif dindex1 == 0 and (pindex2 == 1 or pindex2 == 2):
				print("{:7.3f}".format(ovlp[i-1,j+4]),end="")
				dindex1 = 1
			elif dindex1 == 1 and (pindex2 == 1 or pindex2 == 2):
				print("{:7.3f}".format(ovlp[i-1,j+1]),end="")
				dindex1 = 2
			elif dindex1 == 2 and (pindex2 == 1 or pindex2 == 2):
				print("{:7.3f}".format(ovlp[i-1,j-2]),end="")
				dindex1 = 3
			elif dindex1 == 3 and (pindex2 == 1 or pindex2 == 2):
				print("{:7.3f}".format(ovlp[i-1,j-2]),end="")
				dindex1 = 4
			elif dindex1 == 4 and (pindex2 == 1 or pindex2 == 2):
				print("{:7.3f}".format(ovlp[i-1,j-1]),end="")
				dindex1 = 0
		elif i < ns2 + 3*np2 and j < ns1 + 3*np1 + 5*nd1 + 7*nf1:
			if findex1 == 0 and pindex2 == 0:
				print("{:7.3f}".format(ovlp[i+2,j+6]),end="")
				findex1 = 1
			elif findex1 == 1 and pindex2 == 0:
				print("{:7.3f}".format(ovlp[i+2,j+3]),end="")
				findex1 = 2
			elif findex1 == 2 and pindex2 == 0:
				print("{:7.3f}".format(ovlp[i+2,j]),end="")
				findex1 = 3
			elif findex1 == 3 and pindex2 == 0:
				print("{:7.3f}".format(ovlp[i+2,j-3]),end="")
				findex1 = 4
			elif findex1 == 4 and pindex2 == 0:
				print("{:7.3f}".format(ovlp[i+2,j-3]),end="")
				findex1 = 5
			elif findex1 == 5 and pindex2 == 0:
				print("{:7.3f}".format(ovlp[i+2,j-2]),end="")
				findex1 = 6
			elif findex1 == 6 and pindex2 == 0:
				print("{:7.3f}".format(ovlp[i+2,j-1]),end="")
				findex1 = 0
				pindex2 = 1
			elif findex1 == 0 and (pindex2 == 1 or pindex2 == 2):
				print("{:7.3f}".format(ovlp[i-1,j+6]),end="")
				findex1 = 1
			elif findex1 == 1 and (pindex2 == 1 or pindex2 == 2):
				print("{:7.3f}".format(ovlp[i-1,j+3]),end="")
				findex1 = 2
			elif findex1 == 2 and (pindex2 == 1 or pindex2 == 2):
				print("{:7.3f}".format(ovlp[i-1,j]),end="")
				findex1 = 3
			elif findex1 == 3 and (pindex2 == 1 or pindex2 == 2):
				print("{:7.3f}".format(ovlp[i-1,j-3]),end="")
				findex1 = 4
			elif findex1 == 4 and (pindex2 == 1 or pindex2 == 2):
				print("{:7.3f}".format(ovlp[i-1,j-3]),end="")
				findex1 = 5
			elif findex1 == 5 and (pindex2 == 1 or pindex2 == 2):
				print("{:7.3f}".format(ovlp[i-1,j-2]),end="")
				findex1 = 6
			elif findex1 == 6 and (pindex2 == 1 or pindex2 == 2):
				print("{:7.3f}".format(ovlp[i-1,j-1]),end="")
				findex1 = 0
				if pindex2 == 1:
					pindex2 = 2
				else:
					pindex2 = 0
		elif i < ns2 + 3*np2 + 5*nd2 and j < ns1:
			if dindex2 == 0:
				print("{:7.3f}".format(ovlp[i+4,j]),end="")
			elif dindex2 == 1:
				print("{:7.3f}".format(ovlp[i+1,j]),end="")
			elif dindex2 == 2:
				print("{:7.3f}".format(ovlp[i-2,j]),end="")
			elif dindex2 == 3:
				print("{:7.3f}".format(ovlp[i-2,j]),end="")
			elif dindex2 == 4:
				print("{:7.3f}".format(ovlp[i-1,j]),end="")
		elif i < ns2 + 3*np2 + 5*nd2 and j < ns1 + 3*np1:
			if dindex2 == 0 and pindex1 == 0:
				print("{:7.3f}".format(ovlp[i+4,j+2]),end="")
				pindex1 = 1
			elif dindex2 == 1 and pindex1 == 0:
				print("{:7.3f}".format(ovlp[i+1,j+2]),end="")
				pindex1 = 1
			elif dindex2 == 2 and pindex1 == 0:
				print("{:7.3f}".format(ovlp[i-2,j+2]),end="")
				pindex1 = 1
			elif dindex2 == 3 and pindex1 == 0:
				print("{:7.3f}".format(ovlp[i-2,j+2]),end="")
				pindex1 = 1
			elif dindex2 == 4 and pindex1 == 0:
				print("{:7.3f}".format(ovlp[i-1,j+2]),end="")
				pindex1 = 1
			elif dindex2 == 0 and pindex1 == 1:
				print("{:7.3f}".format(ovlp[i+4,j-1]),end="")
				pindex1 = 2
			elif dindex2 == 1 and pindex1 == 1:
				print("{:7.3f}".format(ovlp[i+1,j-1]),end="")
				pindex1 = 2
			elif dindex2 == 2 and pindex1 == 1:
				print("{:7.3f}".format(ovlp[i-2,j-1]),end="")
				pindex1 = 2
			elif dindex2 == 3 and pindex1 == 1:
				print("{:7.3f}".format(ovlp[i-2,j-1]),end="")
				pindex1 = 2
			elif dindex2 == 4 and pindex1 == 1:
				print("{:7.3f}".format(ovlp[i-1,j-1]),end="")
				pindex1 = 2
			elif dindex2 == 0 and pindex1 == 2:
				print("{:7.3f}".format(ovlp[i+4,j-1]),end="")
				pindex1 = 0
			elif dindex2 == 1 and pindex1 == 2:
				print("{:7.3f}".format(ovlp[i+1,j-1]),end="")
				pindex1 = 0
			elif dindex2 == 2 and pindex1 == 2:
				print("{:7.3f}".format(ovlp[i-2,j-1]),end="")
				pindex1 = 0
			elif dindex2 == 3 and pindex1 == 2:
				print("{:7.3f}".format(ovlp[i-2,j-1]),end="")
				pindex1 = 0
			elif dindex2 == 4 and pindex1 == 2:
				print("{:7.3f}".format(ovlp[i-1,j-1]),end="")
				pindex1 = 0
		elif i < ns2 + 3*np2 + 5*nd2 and j < ns1 + 3*np1 + 5*nd1:
			if dindex1 == 0 and dindex2 == 0:
				print("{:7.3f}".format(ovlp[i+4,j+4]),end="")
				dindex1 = 1
			elif dindex1 == 1 and dindex2 == 0:
				print("{:7.3f}".format(ovlp[i+4,j+1]),end="")
				dindex1 = 2
			elif dindex1 == 2 and dindex2 == 0:
				print("{:7.3f}".format(ovlp[i+4,j-2]),end="")
				dindex1 = 3
			elif dindex1 == 3 and dindex2 == 0:
				print("{:7.3f}".format(ovlp[i+4,j-2]),end="")
				dindex1 = 4
			elif dindex1 == 4 and dindex2 == 0:
				print("{:7.3f}".format(ovlp[i+4,j-1]),end="")
				dindex1 = 0
			elif dindex1 == 0 and dindex2 == 1:
				print("{:7.3f}".format(ovlp[i+1,j+4]),end="")
				dindex1 = 1
			elif dindex1 == 1 and dindex2 == 1:
				print("{:7.3f}".format(ovlp[i+1,j+1]),end="")
				dindex1 = 2
			elif dindex1 == 2 and dindex2 == 1:
				print("{:7.3f}".format(ovlp[i+1,j-2]),end="")
				dindex1 = 3
			elif dindex1 == 3 and dindex2 == 1:
				print("{:7.3f}".format(ovlp[i+1,j-2]),end="")
				dindex1 = 4
			elif dindex1 == 4 and dindex2 == 1:
				print("{:7.3f}".format(ovlp[i+1,j-1]),end="")
				dindex1 = 0
			elif dindex1 == 0 and (dindex2 == 2 or dindex2 == 3):
				print("{:7.3f}".format(ovlp[i-2,j+4]),end="")
				dindex1 = 1
			elif dindex1 == 1 and (dindex2 == 2 or dindex2 == 3):
				print("{:7.3f}".format(ovlp[i-2,j+1]),end="")
				dindex1 = 2
			elif dindex1 == 2 and (dindex2 == 2 or dindex2 == 3):
				print("{:7.3f}".format(ovlp[i-2,j-2]),end="")
				dindex1 = 3
			elif dindex1 == 3 and (dindex2 == 2 or dindex2 == 3):
				print("{:7.3f}".format(ovlp[i-2,j-2]),end="")
				dindex1 = 4
			elif dindex1 == 4 and (dindex2 == 2 or dindex2 == 3):
				print("{:7.3f}".format(ovlp[i-2,j-1]),end="")
				dindex1 = 0
			elif dindex1 == 0 and dindex2 == 4:
				print("{:7.3f}".format(ovlp[i-1,j+4]),end="")
				dindex1 = 1
			elif dindex1 == 1 and dindex2 == 4:
				print("{:7.3f}".format(ovlp[i-1,j+1]),end="")
				dindex1 = 2
			elif dindex1 == 2 and dindex2 == 4:
				print("{:7.3f}".format(ovlp[i-1,j-2]),end="")
				dindex1 = 3
			elif dindex1 == 3 and dindex2 == 4:
				print("{:7.3f}".format(ovlp[i-1,j-2]),end="")
				dindex1 = 4
			elif dindex1 == 4 and dindex2 == 4:
				print("{:7.3f}".format(ovlp[i-1,j-1]),end="")
				dindex1 = 0
		elif i < ns2 + 3*np2 + 5*nd2 and j < ns1 + 3*np1 + 5*nd1 + 7*nf1:
			if findex1 == 0 and dindex2 == 0:
				print("{:7.3f}".format(ovlp[i+4,j+6]),end="")
				findex1 = 1
			elif findex1 == 1 and dindex2 == 0:
				print("{:7.3f}".format(ovlp[i+4,j+3]),end="")
				findex1 = 2
			elif findex1 == 2 and dindex2 == 0:
				print("{:7.3f}".format(ovlp[i+4,j]),end="")
				findex1 = 3
			elif findex1 == 3 and dindex2 == 0:
				print("{:7.3f}".format(ovlp[i+4,j-3]),end="")
				findex1 = 4
			elif findex1 == 4 and dindex2 == 0:
				print("{:7.3f}".format(ovlp[i+4,j-3]),end="")
				findex1 = 5
			elif findex1 == 5 and dindex2 == 0:
				print("{:7.3f}".format(ovlp[i+4,j-2]),end="")
				findex1 = 6
			elif findex1 == 6 and dindex2 == 0:
				print("{:7.3f}".format(ovlp[i+4,j-1]),end="")
				findex1 = 0
				dindex2 = 1
			elif findex1 == 0 and dindex2 == 1:
				print("{:7.3f}".format(ovlp[i+1,j+6]),end="")
				findex1 = 1
			elif findex1 == 1 and dindex2 == 1:
				print("{:7.3f}".format(ovlp[i+1,j+3]),end="")
				findex1 = 2
			elif findex1 == 2 and dindex2 == 1:
				print("{:7.3f}".format(ovlp[i+1,j]),end="")
				findex1 = 3
			elif findex1 == 3 and dindex2 == 1:
				print("{:7.3f}".format(ovlp[i+1,j-3]),end="")
				findex1 = 4
			elif findex1 == 4 and dindex2 == 1:
				print("{:7.3f}".format(ovlp[i+1,j-3]),end="")
				findex1 = 5
			elif findex1 == 5 and dindex2 == 1:
				print("{:7.3f}".format(ovlp[i+1,j-2]),end="")
				findex1 = 6
			elif findex1 == 6 and dindex2 == 1:
				print("{:7.3f}".format(ovlp[i+1,j-1]),end="")
				findex1 = 0
				dindex2 = 2
			elif findex1 == 0 and (dindex2 == 2 or dindex2 == 3):
				print("{:7.3f}".format(ovlp[i-2,j+6]),end="")
				findex1 = 1
			elif findex1 == 1 and (dindex2 == 2 or dindex2 == 3):
				print("{:7.3f}".format(ovlp[i-2,j+3]),end="")
				findex1 = 2
			elif findex1 == 2 and (dindex2 == 2 or dindex2 == 3):
				print("{:7.3f}".format(ovlp[i-2,j]),end="")
				findex1 = 3
			elif findex1 == 3 and (dindex2 == 2 or dindex2 == 3):
				print("{:7.3f}".format(ovlp[i-2,j-3]),end="")
				findex1 = 4
			elif findex1 == 4 and (dindex2 == 2 or dindex2 == 3):
				print("{:7.3f}".format(ovlp[i-2,j-3]),end="")
				findex1 = 5
			elif findex1 == 5 and (dindex2 == 2 or dindex2 == 3):
				print("{:7.3f}".format(ovlp[i-2,j-2]),end="")
				findex1 = 6
			elif findex1 == 6 and (dindex2 == 2 or dindex2 == 3):
				print("{:7.3f}".format(ovlp[i-2,j-1]),end="")
				findex1 = 0
				dindex2 = dindex2 + 1
			elif findex1 == 0 and dindex2 == 4:
				print("{:7.3f}".format(ovlp[i-1,j+6]),end="")
				findex1 = 1
			elif findex1 == 1 and dindex2 == 4:
				print("{:7.3f}".format(ovlp[i-1,j+3]),end="")
				findex1 = 2
			elif findex1 == 2 and dindex2 == 4:
				print("{:7.3f}".format(ovlp[i-1,j]),end="")
				findex1 = 3
			elif findex1 == 3 and dindex2 == 4:
				print("{:7.3f}".format(ovlp[i-1,j-3]),end="")
				findex1 = 4
			elif findex1 == 4 and dindex2 == 4:
				print("{:7.3f}".format(ovlp[i-1,j-3]),end="")
				findex1 = 5
			elif findex1 == 5 and dindex2 == 4:
				print("{:7.3f}".format(ovlp[i-1,j-2]),end="")
				findex1 = 6
			elif findex1 == 6 and dindex2 == 4:
				print("{:7.3f}".format(ovlp[i-1,j-1]),end="")
				findex1 = 0
				dindex2 = 0
		elif i < ns2 + 3*np2 + 5*nd2 + 7*nf2 and j < ns1:
			if findex2 == 0:
				print("{:7.3f}".format(ovlp[i+6,j]),end="")
			elif findex2 == 1:
				print("{:7.3f}".format(ovlp[i+3,j]),end="")
			elif findex2 == 2:
				print("{:7.3f}".format(ovlp[i,j]),end="")
			elif findex2 == 3 or findex2 == 4:
				print("{:7.3f}".format(ovlp[i-3,j]),end="")
			elif findex2 == 5:
				print("{:7.3f}".format(ovlp[i-2,j]),end="")
			elif findex2 == 6:
				print("{:7.3f}".format(ovlp[i-1,j]),end="")
		elif i < ns2 + 3*np2 + 5*nd2 + 7*nf2 and j < ns1 + 3*np1:
			if findex2 == 0 and pindex1 == 0:
				print("{:7.3f}".format(ovlp[i+6,j+2]),end="")
				pindex1 = 1
			elif findex2 == 1 and pindex1 == 0:
				print("{:7.3f}".format(ovlp[i+3,j+2]),end="")
				pindex1 = 1
			elif findex2 == 2 and pindex1 == 0:
				print("{:7.3f}".format(ovlp[i,j+2]),end="")
				pindex1 = 1
			elif (findex2 == 3 or findex2 == 4) and pindex1 == 0:
				print("{:7.3f}".format(ovlp[i-3,j+2]),end="")
				pindex1 = 1
			elif findex2 == 5 and pindex1 == 0:
				print("{:7.3f}".format(ovlp[i-2,j+2]),end="")
				pindex1 = 1
			elif findex2 == 6 and pindex1 == 0:
				print("{:7.3f}".format(ovlp[i-1,j+2]),end="")
				pindex1 = 1
			elif findex2 == 0 and pindex1 == 1:
				print("{:7.3f}".format(ovlp[i+6,j-1]),end="")
				pindex1 = 2
			elif findex2 == 1 and pindex1 == 1:
				print("{:7.3f}".format(ovlp[i+3,j-1]),end="")
				pindex1 = 2
			elif findex2 == 2 and pindex1 == 1:
				print("{:7.3f}".format(ovlp[i,j-1]),end="")
				pindex1 = 2
			elif (findex2 == 3 or findex2 == 4) and pindex1 == 1:
				print("{:7.3f}".format(ovlp[i-3,j-1]),end="")
				pindex1 = 2
			elif findex2 == 5 and pindex1 == 1:
				print("{:7.3f}".format(ovlp[i-2,j-1]),end="")
				pindex1 = 2
			elif findex2 == 6 and pindex1 == 1:
				print("{:7.3f}".format(ovlp[i-1,j-1]),end="")
				pindex1 = 2
			elif findex2 == 0 and pindex1 == 2:
				print("{:7.3f}".format(ovlp[i+6,j-1]),end="")
				pindex1 = 0
			elif findex2 == 1 and pindex1 == 2:
				print("{:7.3f}".format(ovlp[i+3,j-1]),end="")
				pindex1 = 0
			elif findex2 == 2 and pindex1 == 2:
				print("{:7.3f}".format(ovlp[i,j-1]),end="")
				pindex1 = 0
			elif (findex2 == 3 or findex2 == 4) and pindex1 == 2:
				print("{:7.3f}".format(ovlp[i-3,j-1]),end="")
				pindex1 = 0
			elif findex2 == 5 and pindex1 == 2:
				print("{:7.3f}".format(ovlp[i-2,j-1]),end="")
				pindex1 = 0
			elif findex2 == 6 and pindex1 == 2:
				print("{:7.3f}".format(ovlp[i-1,j-1]),end="")
				pindex1 = 0
		elif i < ns2 + 3*np2 + 5*nd2 + 7*nf2 and j < ns1 + 3*np1 + 5*nd1:
			if dindex1 == 0 and findex2 == 0:
				print("{:7.3f}".format(ovlp[i+6,j+4]),end="")
				dindex1 = 1
			elif dindex1 == 1 and findex2 == 0:
				print("{:7.3f}".format(ovlp[i+6,j+1]),end="")
				dindex1 = 2
			elif dindex1 == 2 and findex2 == 0:
				print("{:7.3f}".format(ovlp[i+6,j-2]),end="")
				dindex1 = 3
			elif dindex1 == 3 and findex2 == 0:
				print("{:7.3f}".format(ovlp[i+6,j-2]),end="")
				dindex1 = 4
			elif dindex1 == 4 and findex2 == 0:
				print("{:7.3f}".format(ovlp[i+6,j-1]),end="")
				dindex1 = 0
			elif dindex1 == 0 and findex2 == 1:
				print("{:7.3f}".format(ovlp[i+3,j+4]),end="")
				dindex1 = 1
			elif dindex1 == 1 and findex2 == 1:
				print("{:7.3f}".format(ovlp[i+3,j+1]),end="")
				dindex1 = 2
			elif dindex1 == 2 and findex2 == 1:
				print("{:7.3f}".format(ovlp[i+3,j-2]),end="")
				dindex1 = 3
			elif dindex1 == 3 and findex2 == 1:
				print("{:7.3f}".format(ovlp[i+3,j-2]),end="")
				dindex1 = 4
			elif dindex1 == 4 and findex2 == 1:
				print("{:7.3f}".format(ovlp[i+3,j-1]),end="")
				dindex1 = 0
			elif dindex1 == 0 and findex2 == 2:
				print("{:7.3f}".format(ovlp[i,j+4]),end="")
				dindex1 = 1
			elif dindex1 == 1 and findex2 == 2:
				print("{:7.3f}".format(ovlp[i,j+1]),end="")
				dindex1 = 2
			elif dindex1 == 2 and findex2 == 2:
				print("{:7.3f}".format(ovlp[i,j-2]),end="")
				dindex1 = 3
			elif dindex1 == 3 and findex2 == 2:
				print("{:7.3f}".format(ovlp[i,j-2]),end="")
				dindex1 = 4
			elif dindex1 == 4 and findex2 == 2:
				print("{:7.3f}".format(ovlp[i,j-1]),end="")
				dindex1 = 0
			elif dindex1 == 0 and (findex2 == 3 or findex2 == 4):
				print("{:7.3f}".format(ovlp[i-3,j+4]),end="")
				dindex1 = 1
			elif dindex1 == 1 and (findex2 == 3 or findex2 == 4):
				print("{:7.3f}".format(ovlp[i-3,j+1]),end="")
				dindex1 = 2
			elif dindex1 == 2 and (findex2 == 3 or findex2 == 4):
				print("{:7.3f}".format(ovlp[i-3,j-2]),end="")
				dindex1 = 3
			elif dindex1 == 3 and (findex2 == 3 or findex2 == 4):
				print("{:7.3f}".format(ovlp[i-3,j-2]),end="")
				dindex1 = 4
			elif dindex1 == 4 and (findex2 == 3 or findex2 == 4):
				print("{:7.3f}".format(ovlp[i-3,j-1]),end="")
				dindex1 = 0
			elif dindex1 == 0 and findex2 == 5:
				print("{:7.3f}".format(ovlp[i-2,j+4]),end="")
				dindex1 = 1
			elif dindex1 == 1 and findex2 == 5:
				print("{:7.3f}".format(ovlp[i-2,j+1]),end="")
				dindex1 = 2
			elif dindex1 == 2 and findex2 == 5:
				print("{:7.3f}".format(ovlp[i-2,j-2]),end="")
				dindex1 = 3
			elif dindex1 == 3 and findex2 == 5:
				print("{:7.3f}".format(ovlp[i-2,j-2]),end="")
				dindex1 = 4
			elif dindex1 == 4 and findex2 == 5:
				print("{:7.3f}".format(ovlp[i-2,j-1]),end="")
				dindex1 = 0
			elif dindex1 == 0 and findex2 == 6:
				print("{:7.3f}".format(ovlp[i-1,j+4]),end="")
				dindex1 = 1
			elif dindex1 == 1 and findex2 == 6:
				print("{:7.3f}".format(ovlp[i-1,j+1]),end="")
				dindex1 = 2
			elif dindex1 == 2 and findex2 == 6:
				print("{:7.3f}".format(ovlp[i-1,j-2]),end="")
				dindex1 = 3
			elif dindex1 == 3 and findex2 == 6:
				print("{:7.3f}".format(ovlp[i-1,j-2]),end="")
				dindex1 = 4
			elif dindex1 == 4 and findex2 == 6:
				print("{:7.3f}".format(ovlp[i-1,j-1]),end="")
				dindex1 = 0
		elif i < ns2 + 3*np2 + 5*nd2 + 7*nf2 and j < ns1 + 3*np1 + 5*nd1 + 7*nf1:
			if findex1 == 0 and findex2 == 0:
				print("{:7.3f}".format(ovlp[i+6,j+6]),end="")
				findex1 = 1
			elif findex1 == 1 and findex2 == 0:
				print("{:7.3f}".format(ovlp[i+6,j+3]),end="")
				findex1 = 2
			elif findex1 == 2 and findex2 == 0:
				print("{:7.3f}".format(ovlp[i+6,j]),end="")
				findex1 = 3
			elif (findex1 == 3 or findex1 == 4) and findex2 == 0:
				print("{:7.3f}".format(ovlp[i+6,j-3]),end="")
				findex1 = findex1 + 1
			elif findex1 == 5 and findex2 == 0:
				print("{:7.3f}".format(ovlp[i+6,j-2]),end="")
				findex1 = 6
			elif findex1 == 6 and findex2 == 0:
				print("{:7.3f}".format(ovlp[i+6,j-1]),end="")
				findex2 = 1
				findex1 = 0
			elif findex1 == 0 and findex2 == 1:
				print("{:7.3f}".format(ovlp[i+3,j+6]),end="")
				findex1 = 1
			elif findex1 == 1 and findex2 == 1:
				print("{:7.3f}".format(ovlp[i+3,j+3]),end="")
				findex1 = 2
			elif findex1 == 2 and findex2 == 1:
				print("{:7.3f}".format(ovlp[i+3,j]),end="")
				findex1 = 3
			elif (findex1 == 3 or findex1 == 4) and findex2 == 1:
				print("{:7.3f}".format(ovlp[i+3,j-3]),end="")
				findex1 = findex1 + 1
			elif findex1 == 5 and findex2 == 1:
				print("{:7.3f}".format(ovlp[i+3,j-2]),end="")
				findex1 = 6
			elif findex1 == 6 and findex2 == 1:
				print("{:7.3f}".format(ovlp[i+3,j-1]),end="")
				findex2 = 2
				findex1 = 0
			elif findex1 == 0 and findex2 == 2:
				print("{:7.3f}".format(ovlp[i,j+6]),end="")
				findex1 = 1
			elif findex1 == 1 and findex2 == 2:
				print("{:7.3f}".format(ovlp[i,j+3]),end="")
				findex1 = 2
			elif findex1 == 2 and findex2 == 2:
				print("{:7.3f}".format(ovlp[i,j]),end="")
				findex1 = 3
			elif (findex1 == 3 or findex1 == 4) and findex2 == 2:
				print("{:7.3f}".format(ovlp[i,j-3]),end="")
				findex1 = findex1 + 1
			elif findex1 == 5 and findex2 == 2:
				print("{:7.3f}".format(ovlp[i,j-2]),end="")
				findex1 = 6
			elif findex1 == 6 and findex2 == 2:
				print("{:7.3f}".format(ovlp[i,j-1]),end="")
				findex2 = 3
				findex1 = 0
			elif findex1 == 0 and (findex2 == 3 or findex2 == 4):
				print("{:7.3f}".format(ovlp[i-3,j+6]),end="")
				findex1 = 1
			elif findex1 == 1 and (findex2 == 3 or findex2 == 4):
				print("{:7.3f}".format(ovlp[i-3,j+3]),end="")
				findex1 = 2
			elif findex1 == 2 and (findex2 == 3 or findex2 == 4):
				print("{:7.3f}".format(ovlp[i-3,j]),end="")
				findex1 = 3
			elif (findex1 == 3 or findex1 == 4) and (findex2 == 3 or findex2 == 4):
				print("{:7.3f}".format(ovlp[i-3,j-3]),end="")
				findex1 = findex1 + 1
			elif findex1 == 5 and (findex2 == 3 or findex2 == 4):
				print("{:7.3f}".format(ovlp[i-3,j-2]),end="")
				findex1 = 6
			elif findex1 == 6 and (findex2 == 3 or findex2 == 4):
				print("{:7.3f}".format(ovlp[i-3,j-1]),end="")
				findex2 = findex2 + 1
				findex1 = 0
			elif findex1 == 0 and findex2 == 5:
				print("{:7.3f}".format(ovlp[i-2,j+6]),end="")
				findex1 = 1
			elif findex1 == 1 and findex2 == 5:
				print("{:7.3f}".format(ovlp[i-2,j+3]),end="")
				findex1 = 2
			elif findex1 == 2 and findex2 == 5:
				print("{:7.3f}".format(ovlp[i-2,j]),end="")
				findex1 = 3
			elif (findex1 == 3 or findex1 == 4) and findex2 == 5:
				print("{:7.3f}".format(ovlp[i-2,j-3]),end="")
				findex1 = findex1 + 1
			elif findex1 == 5 and findex2 == 5:
				print("{:7.3f}".format(ovlp[i-2,j-2]),end="")
				findex1 = 6
			elif findex1 == 6 and findex2 == 5:
				print("{:7.3f}".format(ovlp[i-2,j-1]),end="")
				findex2 = 6
				findex1 = 0
			elif findex1 == 0 and findex2 == 6:
				print("{:7.3f}".format(ovlp[i-1,j+6]),end="")
				findex1 = 1
			elif findex1 == 1 and findex2 == 6:
				print("{:7.3f}".format(ovlp[i-1,j+3]),end="")
				findex1 = 2
			elif findex1 == 2 and findex2 == 6:
				print("{:7.3f}".format(ovlp[i-1,j]),end="")
				findex1 = 3
			elif (findex1 == 3 or findex1 == 4) and findex2 == 6:
				print("{:7.3f}".format(ovlp[i-1,j-3]),end="")
				findex1 = findex1 + 1
			elif findex1 == 5 and findex2 == 6:
				print("{:7.3f}".format(ovlp[i-1,j-2]),end="")
				findex1 = 6
			elif findex1 == 6 and findex2 == 6:
				print("{:7.3f}".format(ovlp[i-1,j-1]),end="")
				findex2 = 0
				findex1 = 0
	print("")

