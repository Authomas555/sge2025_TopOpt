#An easy interface to call Advect and Meshdist from Charles Dapogny
from numpy import array,zeros
from subprocess import call,DEVNULL
import os

PATH = '"$HOME/bin/advect"'
TEMP_PATH = '$Home/temp'

def write_sol(path,array,dim_space=2,dim_data=1) -> None:
    '''Solution file writter takes a flat numpy array as input'''
    with open(path,'w') as file:
        file.write('MeshVersionFormatted 1\n')
        file.write('\n')
        file.write('Dimension '+str(dim_space)+'\n')
        file.write('\n')
        file.write('SolAtVertices\n')
        file.write(str(int(len(array)/dim_data))+'\n')
        file.write('1 '+str(dim_data)+'\n')
        file.write('\n')
        array = array.reshape((int(len(array)/dim_data),dim_data))
        for i in range(len(array)):
            for j in range(dim_data): file.write(str(array[i,j])+str(' ')*(j!=dim_data-1))
            file.write('\n')
        file.write('\n')
        file.write('End\n')

def read_sol(path):
    with open(path,'r') as file:
        reader = file.readlines()
        len_data = int(reader[5].rstrip('\n'))
        dim_data = int(reader[6].rstrip('\n').split(' ')[1])
        array = zeros((int(len_data/dim_data),dim_data))
        for i,line in enumerate(reader[8:-2]):
            for j,elt in enumerate(line.rstrip('\n').strip().split(' ')):
                array[i,j] = float(elt)
                
        return array.flatten()

def convert_to_mesh(path,dim_space=2):
    if os.path.exists(path+'.mesh'): os.remove(path+'.mesh')
    with open(path+'.mesh','x') as file_o:
        file_o.write('MeshVersionFormatted 1\n')
        file_o.write('\n')
        file_o.write('Dimension\n')
        file_o.write(str(dim_space)+str('\n'))
        file_o.write('\n')
        with open(path,'r') as file_i:
            reader = file_i.readlines()
            v_start_index = int(reader[0].rstrip('\n'))
            e_start_index = int(reader[v_start_index+1].rstrip('\n'))
            # Write vertex
            file_o.write('Vertices\n')
            file_o.write(reader[0]+str('\n'))
            for i in range(1,v_start_index+1): 
                file_o.write(reader[i].rstrip('\n').strip() + ' 0\n')
            # Write triangle
            file_o.write('\n')
            file_o.write('Triangles\n')
            file_o.write(str(e_start_index)+'\n')
            for i in range(v_start_index+2,v_start_index+e_start_index+2): 
                list_temp = list(filter(lambda x : x != '',reader[i].rstrip('\n').strip().split(' ')))[1:]
                str_temp = list_temp[0] + ' ' + list_temp[1] + ' ' + list_temp[2] + ' 0\n'
                file_o.write(str_temp)
            file_o.write('\n')
            file_o.write('End\n')


class advect:
    
    def __init__(self,phi,vit):
        self.path = PATH
        self.temp_save_dir = TEMP_PATH
        self.windows_path = 'C:/'+TEMP_PATH.split('/mnt/c')[1]
        mesh = phi.space.mesh.ngmesh
        mesh.Export(self.windows_path+'mesh_temp','Neutral Format')
        convert_to_mesh(self.windows_path+'mesh_temp')
        write_sol(self.windows_path+'phi_temp.sol',array(phi.vec),dim_data=1)
        write_sol(self.windows_path+'vit_temp.sol',array(vit.vec),dim_data=2)

    def __advect(self,mesh,sol,vel,sol_out):
        call('wsl '+ 'cd ' + self.temp_save_dir +';' + self.path+' '+mesh+' '+"-c "+sol+" "+"-s "+vel+" "+"-o "+sol_out + ' -v')

    def __redist(self,mesh,sol):
        call('wsl '+ 'cd ' + self.temp_save_dir+';'+'cp '+sol+' '+'mesh_temp.sol')
        call('wsl '+ 'cd ' + self.temp_save_dir +';' + 'mshdist '+ mesh)

    def Do(self,redist=False):
        self.__advect('mesh_temp.mesh','phi_temp.sol','vit_temp.sol','phi_out.sol')
        if redist : self.__redist('mesh_temp.mesh','phi_out.sol')
        return read_sol(self.windows_path+'phi_out.sol')

class redist:
    
    def __init__(self,phi):
        self.path = PATH
        self.temp_save_dir = TEMP_PATH
        self.windows_path = 'C:/'+TEMP_PATH.split('/mnt/c')[1]
        mesh = phi.space.mesh.ngmesh
        mesh.Export(self.windows_path+'mesh_temp','Neutral Format')
        convert_to_mesh(self.windows_path+'mesh_temp')
        write_sol(self.windows_path+'phi_temp.sol',array(phi.vec),dim_data=1)

    def __redist(self,mesh,sol):
        call('wsl '+ 'cd ' + self.temp_save_dir+';'+'cp '+sol+' '+'mesh_temp.sol',stderr=DEVNULL, stdout=DEVNULL)
        call('wsl '+ 'cd ' + self.temp_save_dir +';' + 'mshdist '+ mesh,stderr=DEVNULL, stdout=DEVNULL)

    def Do(self):
        self.__redist('mesh_temp.mesh','phi_temp.sol')
        return read_sol(self.windows_path+'mesh_temp.sol')

