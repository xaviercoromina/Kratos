__all__ = ["MPIUnion"]

# Core imports
import KratosMultiphysics

# STD imports
import typing


def GetBroadcastOperation(item, data_communicator: KratosMultiphysics.DataCommunicator) -> typing.Callable:
    if isinstance(item, (int, float, str)):
        return data_communicator.Broadcast
    elif hasattr(item, "__iter__"):
        # Don't allow heterogeneous containers
        types = set([type(component) for component in item])
        if 1 < len(types):
            raise TypeError("Heterogeneous container with types: {types}")

        component_type = types.pop()
        if issubclass(component_type, int):
            return data_communicator.BroadcastInts
        elif issubclass(component_type, float):
            return data_communicator.BroadcastDoubles
        elif issubclass(component_type, str):
            return data_communicator.BroadcastStrings
        else:
            raise TypeError(f"Container of unsupported items: {component_type}")
    else:
        raise TypeError(f"Unsupported type for Broadcast operation: {type(item)}")


def MPIUnion(container: set, data_communicator: KratosMultiphysics.DataCommunicator) -> set:
    union = set()
    if data_communicator.IsDistributed():
        # Send containers from all ranks to the main one
        for rank in range(data_communicator.Size()):
            print(f"Broadcast on {rank} from {data_communicator.Rank()}")
            union |= set(GetBroadcastOperation(container, data_communicator)(container, rank))
    else:
        union = set(container)
    return union
