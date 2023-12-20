with open('def2tzvp.out') as inFile:
	lines = inFile.readlines()

index = 0
while lines[index].find("Basis Dimension") == -1:
	index = index +1

dims = int(lines[index].split()[4])

print(dims)

index = 0
while lines[index].find("groups of distinct atoms") == -1:
	index = index + 1

num_groups = int(lines[index].split()[2])

print(num_groups)
