from KratosMultiphysics import PoromechanicsApplication
import KratosMultiphysics
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
import h5py

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PoromechanicsL2NormErrorProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"

class PoromechanicsL2NormErrorProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters( """
            {
                "model_part_name": "PLEASE_CHOOSE_MODEL_PART_NAME",
                "reference_solution_file_name": "reference_solution.hdf5",
                "write_reference_solution": true
            }  """ )
        settings.ValidateAndAssignDefaults(default_settings)

        self.write_reference_solution = settings["write_reference_solution"].GetBool()

        self.model_part = Model[settings["model_part_name"].GetString()]

        file_name = settings["reference_solution_file_name"].GetString()
        self.file = 0
        if self.write_reference_solution:
            self.file = h5py.File(file_name, 'w')
        else:
            self.file = h5py.File(file_name, 'r')

        self.params = KratosMultiphysics.Parameters("""{}""")
        self.params.AddValue("model_part_name",settings["model_part_name"])
        self.process = KratosPoro.PoromechanicsL2NormErrorProcess(self.model_part, self.params)


    def CreateGroup(self, file_or_group, name, overwrite_previous = True):
        if name in file_or_group:
            if overwrite_previous:
                file_or_group['/'].__delitem__(name)
            else:
                return file_or_group['/' + name]
        return file_or_group.create_group(name)

    def DeleteDataSet(self, file_or_group, dset_name):
        if dset_name in file_or_group:
            file_or_group.__delitem__(dset_name)

    def WriteDataToFile(self, file_or_group, names, data, dtype):
        self.sub_group = self.CreateGroup(file_or_group, self.group_name)
        for name, datum in zip(names, data):
            self.DeleteDataSet(file_or_group, name)
        for name, datum in zip(names, data):
            self.sub_group.create_dataset(name = name, data = datum, dtype = dtype)


    def ExecuteInitialize(self):

        if not self.write_reference_solution:
            self.process.ExecuteInitialize()


    def ExecuteAfterOutputStep(self):

        # TODO
        # time_string = str(round(self.model_part.ProcessInfo[KratosMultiphysics.TIME],12))
        time_string = str(round(self.model_part.ProcessInfo[KratosMultiphysics.TIME],3))
        self.group_name = time_string

        if self.write_reference_solution:
            self.displ_x_list = []
            self.displ_y_list = []
            self.displ_z_list = []
            self.fluid_pressure_list = []
            self.node_id_list = []

            for node in self.model_part.Nodes:
                node_id = node.Id
                displ_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                displ_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                displ_z = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z)
                fluid_pressure = node.GetSolutionStepValue(KratosMultiphysics.WATER_PRESSURE)

                self.node_id_list.append(node_id)
                self.displ_x_list.append(displ_x)
                self.displ_y_list.append(displ_y)
                self.displ_z_list.append(displ_z)
                self.fluid_pressure_list.append(fluid_pressure)

            self.WriteDataToFile(file_or_group = self.file,
                                names = ['NODE', 'REF_DISPL_X', 'REF_DISPL_Y', 'REF_DISPL_Z','REF_FLUID_PRESSURE'],
                                data = [self.node_id_list, self.displ_x_list, self.displ_y_list, self.displ_z_list,self.fluid_pressure_list],
                                dtype = 'float32')
        else:
            e = False
            node = self.group_name + '/REF_DISPL_X'
            if node in self.file.keys():
                # self.file[node]
                e = True
            if e:
                i = 0

                ref_displ_x = self.file[self.group_name + '/REF_DISPL_X'][:,]
                ref_displ_y = self.file[self.group_name + '/REF_DISPL_Y'][:,]
                ref_displ_z = self.file[self.group_name + '/REF_DISPL_Z'][:,]
                ref_fluid_pressure = self.file[self.group_name + '/REF_FLUID_PRESSURE'][:,]

                for node in self.model_part.Nodes:
                    node.SetSolutionStepValue(KratosPoro.REFERENCE_DISPLACEMENT_X,ref_displ_x[i])
                    node.SetSolutionStepValue(KratosPoro.REFERENCE_DISPLACEMENT_Y,ref_displ_y[i])
                    node.SetSolutionStepValue(KratosPoro.REFERENCE_DISPLACEMENT_Z,ref_displ_z[i])
                    node.SetSolutionStepValue(KratosPoro.REFERENCE_FLUID_PRESSURE,ref_fluid_pressure[i])

                    i += 1
                
                self.process.ExecuteAfterOutputStep()