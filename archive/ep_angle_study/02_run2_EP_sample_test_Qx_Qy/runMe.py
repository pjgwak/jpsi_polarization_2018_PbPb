import subprocess as sp

# MC
name_label = 'Run2MC'
nevt = 1000000
isMC = True
command = ['root', '-l', '-b', '-q', f'epTest.C("{name_label}", {nevt}, {int(isMC)})' ]
sp.run(command)


# Data - Cent
name_label = 'Run2DataCent'
nevt = 1000000
isMC = False
isPeri = False
command = ['root', '-l', '-b', '-q', f'epTest.C("{name_label}", {nevt}, {int(isMC)}, {int(isPeri)})' ]
sp.run(command)

# Data - Peri
name_label = 'Run2DataPeri'
nevt = 1000000
isMC = False
isPeri = True
command = ['root', '-l', '-b', '-q', f'epTest.C("{name_label}", {nevt}, {int(isMC)}, {int(isPeri)})' ]
sp.run(command)