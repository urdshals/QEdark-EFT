import os


def change_line_in_file(file, line_begining, replacement):
    a_file = open(file, "r")
    list_of_lines = a_file.readlines()
    for i in range(len(list_of_lines)):
        if list_of_lines[i].startswith(line_begining):
            list_of_lines[i] = replacement
    a_file = open(file, "w")
    a_file.writelines(list_of_lines)
    a_file.close()


def remove_lines_from_file(file, line_begining):
    a_file = open(file, "r")
    list_of_lines = a_file.readlines()
    for i in range(len(list_of_lines)):
        if list_of_lines[i].startswith(line_begining):
            list_of_lines[i] = ''
    a_file = open(file, "w")
    a_file.writelines(list_of_lines)
    a_file.close()


for tubes in ['False']:
    for coupling_nr in [16]:
        for interaction in ['l']:
            os.system('cp all_rate_computation.py all_rate_computation' + str(coupling_nr) + '_tubes_' + tubes + '_interaction' + interaction + '.py')
            os.system('cp run.sh run' + str(coupling_nr) + '_tubes_' + tubes + '_interaction' + interaction + '.sh')
            change_line_in_file('all_rate_computation' + str(coupling_nr) + '_tubes_' + tubes + '_interaction' + interaction + '.py', '    tubes', '    tubes = ' + tubes +'\n')
            change_line_in_file('all_rate_computation' + str(coupling_nr) + '_tubes_' + tubes + '_interaction' + interaction + '.py', '    runname', '    runname = ' +'"'+"all_rates" + str(coupling_nr) + "_tubes" + str(tubes) + "_interaction" + interaction +'"' +'\n')
            change_line_in_file('all_rate_computation' + str(coupling_nr) + '_tubes_' + tubes + '_interaction' + interaction + '.py', '            for interaction', '            for interaction in ["' + interaction + '"]:\n')
            change_line_in_file('all_rate_computation' + str(coupling_nr) + '_tubes_' + tubes + '_interaction' + interaction + '.py', '    for coupling_nr in', '    for coupling_nr in [' + str(coupling_nr) + ']:\n')
            change_line_in_file('run' + str(coupling_nr) + '_tubes_' + tubes + '_interaction' + interaction + '.sh', '#SBATCH -J', '#SBATCH -J '+str(coupling_nr)+interaction+"\n")
            change_line_in_file('run' + str(coupling_nr) + '_tubes_' + tubes + '_interaction' + interaction + '.sh', 'python3', 'python3 ' + 'all_rate_computation' + str(coupling_nr) + '_tubes_' + tubes + '_interaction' + interaction + '.py')
            os.system('sbatch run' + str(coupling_nr) + '_tubes_' + tubes + '_interaction' + interaction + '.sh')

