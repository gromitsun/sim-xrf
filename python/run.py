import os
arglist = ['wpw90_%s' % i for i in range(1,5)]
# arglist = ['as3']
os.chdir('./core/')
os.system('python sim.py %s' % ' '.join(arglist))