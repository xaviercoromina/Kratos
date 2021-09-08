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
                "output_control_type": "step",
                "output_interval": 1.0,
                "write_mode": true
            }  """ )
        settings.ValidateAndAssignDefaults(default_settings)

        self.write_mode = settings["write_mode"].GetBool()

        self.model_part = Model[settings["model_part_name"].GetString()]

        output_control_type = settings["output_control_type"].GetString()
        if output_control_type == "step":
            self.output_control_is_time = False
        else:
            self.output_control_is_time = True
        
        self.output_interval = settings["output_interval"].GetDouble()

        self.step_count = 0
        self.next_output = 0.0

        self.file = 0
        if self.write_mode:
            self.file = h5py.File("reference_displacement.hdf5", 'w')
        else:
            self.file = h5py.File("reference_displacement.hdf5", 'r')

        self.params = KratosMultiphysics.Parameters("""{}""")
        self.params.AddValue("model_part_name",settings["model_part_name"])
        self.process = KratosPoro.PoromechanicsL2NormErrorProcess(self.model_part, self.params)

    def IsOutputStep(self):
        if self.output_control_is_time:
            time = self.__get_pretty_time(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
            return (time >= self.__get_pretty_time(self.next_output))
        else:
            return ( self.step_count >= self.next_output )

    def __get_pretty_time(self,time):
        pretty_time = "{0:.12g}".format(time)
        pretty_time = float(pretty_time)
        return pretty_time

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

    # def ExecuteInitialize(self):

    #     self.process.ExecuteInitialize()

    # def ExecuteInitializeSolutionStep(self):

    #     self.process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):

        # Schedule next output
        time = self.__get_pretty_time(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
        if self.output_interval > 0.0: # Note: if == 0, we'll just always print
            if self.output_control_is_time:
                while self.__get_pretty_time(self.next_output) <= time:
                    self.next_output += self.output_interval
            else:
                while self.next_output <= self.step_count:
                    self.next_output += self.output_interval

        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        time_str = str(round(time,3))
        self.group_name = time_str

        if self.write_mode:
            self.displ_x_list = []
            self.displ_y_list = []
            self.displ_z_list = []
            self.node_id_list = []

            for node in self.model_part.Nodes:
                node_id = node.Id
                displ_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                displ_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                displ_z = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z)

                self.node_id_list.append(node_id)
                self.displ_x_list.append(displ_x)
                self.displ_y_list.append(displ_y)
                self.displ_z_list.append(displ_z)

            self.WriteDataToFile(file_or_group = self.file,
                                names = ['NODE', 'REF_DISPL_X', 'REF_DISPL_Y', 'REF_DISPL_Z'],
                                data = [self.node_id_list, self.displ_x_list, self.displ_y_list, self.displ_z_list],
                                dtype = 'float32')
        else:
            i = 0
            ref_displ_x = self.file[self.group_name + '/REF_DISPL_X'][:,]
            ref_displ_y = self.file[self.group_name + '/REF_DISPL_Y'][:,]
            ref_displ_z = self.file[self.group_name + '/REF_DISPL_Z'][:,]

            for node in self.model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.REFERENCE_DISPLACEMENT_X,ref_displ_x[i])
                node.SetSolutionStepValue(KratosMultiphysics.REFERENCE_DISPLACEMENT_Y,ref_displ_y[i])
                node.SetSolutionStepValue(KratosMultiphysics.REFERENCE_DISPLACEMENT_Z,ref_displ_z[i])

                i += 1
            
            self.process.ExecuteFinalizeSolutionStep()
            # TODO: could we do it directly from python ?
            
