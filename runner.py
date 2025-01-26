import subprocess, os, sys
import datetime
from pathlib import Path

from ExperimentSet import num_qubit_exp, error_rate_exp, trap_size_exp, sched_exp

_filepath = os.path.abspath(__file__)
_dirname = os.path.dirname(_filepath)
_executable = os.path.join(_dirname, "target", "release")
if os.path.isfile(os.path.join(_executable, "qecc.exe")):
    _executable = os.path.join(_executable, "qecc.exe")
elif os.path.isfile(os.path.join(_executable, "qecc")):
    _executable = os.path.join(_executable, "qecc")
else:
    raise Exception(f'No executable exists')

NUM_ITER = 1

def num_qubit_sim():
    for exp in num_qubit_exp:
        for num_qubit, approx_factor in exp.num_qubits:
            filename_qec = f'{exp.algorithm}({num_qubit})_af({approx_factor})_qec'
            filename_noqec = f'{exp.algorithm}({num_qubit})_af({approx_factor})_noqec'

            ### qec
            my_env = os.environ.copy()
            my_env['TRAP_SIZE'] = str(exp.trap_size)
            my_env['INDIV_EXP'] = 'True'
            my_env['NUM_ITER'] = str(NUM_ITER)
            my_env['QEC'] = 'True'
            with open(f'circuit/{filename_qec}.txt', 'r') as f1:
                with open(f'data/num_qubit/{filename_qec}.txt', 'w') as f2:
                    print("%s Start num_qubit %s(%d)_af(%d) with qec"%\
                          (datetime.datetime.today().strftime("[%H:%M:%S]"), exp.algorithm, num_qubit, approx_factor))
                    subprocess.run(_executable, stdin=f1, stdout=f2, env=my_env)
                    print("%s Finish num_qubit %s(%d)_af(%d) with qec"%\
                          (datetime.datetime.today().strftime("[%H:%M:%S]"), exp.algorithm, num_qubit, approx_factor))

            ### no qec
            my_env = os.environ.copy()
            my_env['TRAP_SIZE'] = str(exp.trap_size)
            my_env['INDIV_EXP'] = 'True'
            my_env['NUM_ITER'] = str(NUM_ITER)
            my_env['NO_QEC'] = 'True'
            with open(f'circuit/{filename_noqec}.txt', 'r') as f1:
                with open(f'data/num_qubit/{filename_noqec}.txt', 'w') as f2:
                    print("%s Start num_qubit %s(%d)_af(%d) without qec"%\
                          (datetime.datetime.today().strftime("[%H:%M:%S]"), exp.algorithm, num_qubit, approx_factor))
                    subprocess.run(_executable, stdin=f1, stdout=f2, env=my_env)
                    print("%s Finish num_qubit %s(%d)_af(%d) without qec"%\
                          (datetime.datetime.today().strftime("[%H:%M:%S]"), exp.algorithm, num_qubit, approx_factor))

def error_rate_sim():
    for exp in error_rate_exp:
        num_qubit, approx_factor = exp.num_qubits[0]
        filename_qec = f'{exp.algorithm}({num_qubit})_af({approx_factor})_qec'
        filename_noqec = f'{exp.algorithm}({num_qubit})_af({approx_factor})_noqec'

        ### qec
        my_env = os.environ.copy()
        my_env['TRAP_SIZE'] = str(exp.trap_size)
        my_env['RATE_EXP'] = 'True'
        my_env['NUM_ITER'] = str(NUM_ITER)
        my_env['QEC'] = 'True'
        with open(f'circuit/{filename_qec}.txt', 'r') as f1:
            with open(f'data/error_rate/{filename_qec}_er.txt', 'w') as f2:
                print("%s Start err_rate %s(%d)_af(%d) with qec"%\
                      (datetime.datetime.today().strftime("[%H:%M:%S]"),exp.algorithm, num_qubit, approx_factor))
                subprocess.run(_executable, stdin=f1, stdout=f2, env=my_env)
                print("%s Finish err_rate %s(%d)_af(%d) with qec"%\
                      (datetime.datetime.today().strftime("[%H:%M:%S]"), exp.algorithm, num_qubit, approx_factor))

        ### no qec
        my_env = os.environ.copy()
        my_env['TRAP_SIZE'] = str(exp.trap_size)
        my_env['RATE_EXP'] = 'True'
        my_env['NUM_ITER'] = str(NUM_ITER)
        my_env['NO_QEC'] = 'True'
        with open(f'circuit/{filename_noqec}.txt', 'r') as f1:
            with open(f'data/error_rate/{filename_noqec}_er.txt', 'w') as f2:
                print("%s Start err_rate %s(%d)_af(%d) without qec"%\
                      (datetime.datetime.today().strftime("[%H:%M:%S]"), exp.algorithm, num_qubit, approx_factor))
                subprocess.run(_executable, stdin=f1, stdout=f2, env=my_env)
                print("%s Finish err_rate %s(%d)_af(%d) without qec"%\
                      (datetime.datetime.today().strftime("[%H:%M:%S]"), exp.algorithm, num_qubit, approx_factor))

def trap_size_sim():
    for exp in trap_size_exp:
        num_qubit, approx_factor = exp.num_qubits[0]
        filename_qec = f'{exp.algorithm}({num_qubit})_af({approx_factor})_qec'
        filename_noqec = f'{exp.algorithm}({num_qubit})_af({approx_factor})_noqec'

        ### qec
        my_env = os.environ.copy()
        my_env['SIZE_EXP'] = 'True'
        my_env['NUM_ITER'] = str(NUM_ITER)
        my_env['QEC'] = 'True'
        with open(f'circuit/{filename_qec}.txt', 'r') as f1:
            with open(f'data/trap_size/{filename_qec}_ss.txt', 'w') as f2:
                print("%s Start trap_size %s(%d)_af(%d) with qec"%\
                      (datetime.datetime.today().strftime("[%H:%M:%S]"), exp.algorithm, num_qubit, approx_factor))
                subprocess.run(_executable, stdin=f1, stdout=f2, env=my_env)
                print("%s Finish trap_size %s(%d)_af(%d) with qec"%\
                      (datetime.datetime.today().strftime("[%H:%M:%S]"), exp.algorithm, num_qubit, approx_factor))

        ### no qec
        my_env = os.environ.copy()
        my_env['SIZE_EXP'] = 'True'
        my_env['NUM_ITER'] = str(NUM_ITER)
        my_env['NO_QEC'] = 'True'
        with open(f'circuit/{filename_noqec}.txt', 'r') as f1:
            with open(f'data/trap_size/{filename_noqec}_ss.txt', 'w') as f2:
                print("%s Start trap_size %s(%d)_af(%d) without qec"%\
                      (datetime.datetime.today().strftime("[%H:%M:%S]"), exp.algorithm, num_qubit, approx_factor))
                subprocess.run(_executable, stdin=f1, stdout=f2, env=my_env)
                print("%s Finish trap_size %s(%d)_af(%d) without qec"%\
                      (datetime.datetime.today().strftime("[%H:%M:%S]"), exp.algorithm, num_qubit, approx_factor))

def sched_sim():
    for exp in sched_exp:
        num_qubit, approx_factor = exp.num_qubits[0]
        trap_size = str(exp.trap_size)
        filename_qec = f'{exp.algorithm}({num_qubit})_af({approx_factor})_qec'

        ### qec
        my_env = os.environ.copy()
        my_env['TRAP_SIZE'] = str(exp.trap_size)
        my_env['SCHED_EXP'] = 'True'
        with open(f'circuit/{filename_qec}.txt', 'r') as f1:
            with open(f'data/sched/{filename_qec}_sched_ss{trap_size}.txt', 'w') as f2:
                print("%s Start num_qubit %s(%d)_af(%d) with qec"%\
                        (datetime.datetime.today().strftime("[%H:%M:%S]"), exp.algorithm, num_qubit, approx_factor))
                subprocess.run(_executable, stdin=f1, stdout=f2, env=my_env)
                print("%s Finish num_qubit %s(%d)_af(%d) with qec"%\
                        (datetime.datetime.today().strftime("[%H:%M:%S]"), exp.algorithm, num_qubit, approx_factor))

if __name__ == '__main__':
    exp_name = sys.argv[1]
    if exp_name=='er':
        error_rate_sim()
    elif exp_name=='ts':
        trap_size_sim()
    elif exp_name=='qn':
        num_qubit_sim()
    elif exp_name=="sched":
        sched_sim()
    else:
        raise Exception(f'Invalid command %s - choose among er, ss, qn'%exp_name)