# Core imports
import KratosMultiphysics

# STD imports
import pkgutil
import importlib


class ImportUtility:
    """@brief Utility class for importing modules/objects at runtime."""

    def __init__(self, parameters: KratosMultiphysics.Parameters):
        """@brief Construct an ImportUtility instance to help with parsing @ref Parameters objects.
           @details The only required key in the input @ref Parameters is 'import_name', which is the
                    name of the object to be imported.
           @details An optional key 'import_module' can be specified
                    to narrow down the search (the default value 'KratosMultiphysics' leads to
                    every submodule/application getting imported).
           @details The optional key "max_recursion_depth" limits how many levels of submodules
                    are searched. The search defaults to unlimited recursion.
           @details All other keys are ignored.
           @details For example, the following configuration imports this class:
                    {
                        "import_name" : "ImportUtility"                             <== name of the object to import
                        "import_module" : "KratosMultiphysics.HDF5Application"      <== a guess which module the desired object might be part of
                        "import_recursion_depth" : 1,                               <== 1 level of submodules is searched
                        "some_other_setting" : "some_value"                         <== ignored, no exception raised
                    }
        """
        if not parameters.Has("import_name"):
            raise ValueError(f"Expecting parameters containing 'import_name', got {parameters}")
        self.__object_name = parameters["import_name"]

        if parameters.Has("import_module"):
            self.__module_name = parameters["import_module"]
        else:
            self.__module_name = "KratosMutliphysics"

        if parameters.Has("import_recursion_depth"):
            self.__max_recursion_depth = parameters["import_recursion_depth"].GetInt()
        else:
            self.__max_recursion_depth = None

    def __call__(self) -> type:
        """@brief Get the object specified in the input parameters."""
        return self.GetObject(self.__object_name, module_name = self.__module_name, max_recursion_depth = self.__max_recursion_depth)

    @staticmethod
    def GetObject(object_name: str,
                  module_name: str = "KratosMultiphysics",
                  max_recursion_depth: int = None) -> type:
        """@brief Search a module for an object by name.
           @param object_name: name of the object to import (without module scope).
           @param module_name: full module/submodule path to search.
           @param max_recursion_depth: max depth of nested submodules to search (default is unlimited).
           @details All submodules of the provided module are searched for the requested
                    object. An exception is thrown if multiple objects with the same
                    name are found.
           @warning The module and all encountered submodules get imported during the search.
        """
        objects = []
        ImportUtility._CollectObjects(module_name, object_name, objects, max_recursion_depth = max_recursion_depth)

        # Get rid of duplicates
        ids = set()
        objects = [ids.add(id(obj["object"])) or obj for obj in objects if not (id(obj["object"]) in ids)]

        # Check for objects with identical names
        if not objects:
            raise RuntimeError(f"No object named '{object_name}' found in module '{module_name}'")
        elif 1 < len(objects):
            new_line = "\n"
            for o in objects:
                print(id(o["object"]))
            raise RuntimeError(f"Multiple definitions of '{object_name}' found in module '{module_name}':{new_line}{new_line.join(hit['module_name'] for hit in objects)}")

        return objects[0]["object"]

    @staticmethod
    def _CollectObjects(module_name: str,
                        object_name: str,
                        objects: list,
                        max_recursion_depth: int,
                        depth: int = 0) -> None:
        """@brief Recursively collect all objects from a module matching the provided name.
           @param module_name: full module path to search in.
           @param object_name: name of the object to search for.
           @param max_recursion_depth: max depth of nested submodules to search (default is unlimited).
           @details The found objects get appended to the passed list.
        """
        module = ImportUtility._ImportModule(module_name)
        for item in pkgutil.iter_modules(module.__path__):
            item_name = f"{module_name}.{item.name}"
            if item.ispkg and (max_recursion_depth == None or max_recursion_depth < depth):
                ImportUtility._CollectObjects(item_name, object_name, objects, max_recursion_depth, depth = depth + 1)
            else:
                item = ImportUtility._ImportModule(item_name)
                for member_name in dir(item):
                    if member_name == object_name:
                        objects.append({
                            "object" : getattr(item, member_name),
                            "module_name" : item_name
                        })

    @staticmethod
    def _ImportModule(module_name: str) -> "module":
        try:
            return importlib.import_module(module_name)
        except Exception as exception:
            KratosMultiphysics.Logger.PrintWarning(f"[ImportUtility] Importing module '{module_name}' failed with message: {exception}")
