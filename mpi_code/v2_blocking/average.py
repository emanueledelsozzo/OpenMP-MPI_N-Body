import sys

filename = sys.argv[1]

initGrid = "InitGrid"
initGridVal = []
genField = "GenerationField"
genFieldVal = []
partGen = "ParticleGeneration"
partGenVal = []
sysEv = "SystemEvolution"
sysEvVal = []
node_compt = "node_computation"
node_comptVal = []
node_comm = "node_communication"
node_commVal = []
init_mpi = "InitMPI"
init_mpiVal = []

f = open(filename, "r")
for i in f.readlines():
	tok = i.split(" ")
	if initGrid == tok[0]:
		initGridVal.append(float(tok[-2]))
	elif genField == tok[0]:
		genFieldVal.append(float(tok[-2]))
	elif partGen == tok[0]:
		partGenVal.append(float(tok[-2]))
	elif sysEv == tok[0]:
		sysEvVal.append(float(tok[-2]))
	elif node_compt == tok[0]:
		node_comptVal.append(float(tok[-2]))
	elif node_comm == tok[0]:
                node_commVal.append(float(tok[-2]))
	elif init_mpi == tok[0]:
                init_mpiVal.append(float(tok[-2]))

print "initGrid list = %s" % (str(initGridVal))
print "genField list = %s" % (str(genFieldVal))
print "partGen list  = %s" % (str(partGenVal))
print "sysEv list    = %s" % (str(sysEvVal))
print "node compt    = %s" % (str(node_comptVal))
print "node comm     = %s" % (str(node_commVal))
print "init MPI      = %s" % (str(init_mpiVal))

print "initGrid   = %.10f" % (sum(initGridVal)/len(initGridVal))
print "genField   = %.10f" % (sum(genFieldVal)/len(genFieldVal))
print "partGen    = %.10f" % (sum(partGenVal)/len(partGenVal))
print "sysEv      = %.10f" % (sum(sysEvVal)/len(sysEvVal))
print "node compt = %.10f" % (sum(node_comptVal)/len(node_comptVal))
print "node comm  = %.10f" % (sum(node_commVal)/len(node_commVal))
print "init MPI   = %.10f" % (sum(init_mpiVal)/len(init_mpiVal))
