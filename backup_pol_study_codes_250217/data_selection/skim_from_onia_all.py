import subprocess


date = '990000'

# Check ../cms_headers/cutsAndBin.h
kTrigJpsi = 12

# start time
print(kTrigJpsi)

print('Start Data DM')
command = f"root -l -b -q onia_to_skim_jpsi.C'(100, false, 1, {kTrigJpsi}, 0, 1, \"{date}\")'"
subprocess.run(command, shell=True)
print('End Data DM')

#print('Start Data Peri')
#print('End Data Peri')

#print('Start MC Signal')
#print('End MC Signal')


#print('Start MC BtoJpsi')
#print('Start MC BtoJpsi')



# total time