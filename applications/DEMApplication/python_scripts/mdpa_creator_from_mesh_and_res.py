import KratosMultiphysics
from KratosMultiphysics.DEMApplication import mesh_creator_sphere_2D
import os

problem_name = "ctw10_20"
main_path = os.getcwd()
spheres_mp_filename_post = problem_name + 'DEM_Post'
desired_time_to_print = '0.001026'
file_msh = problem_name + '_' + desired_time_to_print + '.post.msh'
file_res = problem_name + '_' + desired_time_to_print + '.post.res'
post_path = main_path

mesh_creator_sphere_2D.WriteSphereMdpaFromResults(problem_name + 'DEM', main_path, spheres_mp_filename_post, file_msh, file_res, post_path)
