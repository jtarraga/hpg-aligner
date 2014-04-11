import subprocess
import threading
import sys

num_machines = 1
if len(sys.argv) > 1:
    num_machines = int(sys.argv[1])

if (num_machines <= 0):
    num_machines = 1;

"""
dataset = "/home/jtarraga/datasets/chr21_100nt_r0.01.bwa.read1.fastq"
"""
dataset = "/data/hpg/fastq/simulated/human/DNA/GRCh37_ens73/4M_100nt_r0.01.bwa.read1.fastq"

hpg_aligner = "/home/jtarraga/appl/bioinfo-c/hpg-aligner/bin/hpg-aligner dna -f " + dataset + " -i /tmp/sa-index-chr20-21-22/ --report-n-best=1 --cpu-threads 12 -o /home/jtarraga/appl/bioinfo-c/hpg-aligner/hpg_aligner_out" 

threads = []
machines = ["clinic10", "clinic11", "clinic12"]
machine_id = 0

def thread_function(machine, pid, np):
    conn = "jtarraga@" + machine
    cmdline = hpg_aligner + " --id " + str(pid) + " --np " + str(np)
    subprocess.call(["ssh", conn, cmdline])
    
print "Creating " + str(num_machines) +  " threads"
for i in range (0, num_machines):
    machine = machines[i]
    thread = threading.Thread(target=thread_function, args=(machine, machine_id, num_machines, ))
    threads += [thread]
    machine_id += 1

print "Launching threads"
for thread in threads:
    thread.start()

print "Waiting for threads"
for thread in threads:
    thread.join()

print "Done."

