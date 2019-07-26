import subprocess
import time
MAX_JOBS = 75

# Get number of jobs in queue
num_jobs = subprocess.check_output('squeue -u USERNAME | wc -l',shell=True)
num_jobs = int(num_jobs) - 1

# Wait until opening in queue    
while (num_jobs >= MAX_JOBS): 
	time.sleep(60)
	num_jobs = subprocess.check_output('squeue -u USERNAME | wc -l',shell=True)
	num_jobs = int(num_jobs) - 1