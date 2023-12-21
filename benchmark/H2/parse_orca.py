with open('def2tzvp.out') as inFile:
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
	bs1_dim = dims/2
	bs2_dim = dims/2
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
		print(str2[:sloc2])
	ploc2 = str2.find("p")
	if ploc2 != -1:
		print(str2[sloc2+1:ploc2])
	dloc2 = str2.find("d")
	if dloc2 != -1:
		print(str2[ploc2+1:dloc2])
	floc2 = str2.find("f")
	if floc2 != -1:
		print(str2[dloc2+1:floc2])
	gloc2 = str2.find("g")
	if gloc2 != -1:
		print(str2[floc2+1:gloc2])
	hloc2 = str2.find("h")
	if hloc2 != -1:
		print(str2[gloc2+1:hloc2])
	iloc2 = str2.find("i")
	if iloc2 != -1:
		print(str2[hloc2+1:iloc2])


print(ns1+3*np1)
