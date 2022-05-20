---
title: Source Documentation
keywords: documentation, doxygen
tags: [Source_Documentation.md]
sidebar: kratos_for_developers
summary:
---
## Intro

This page is meant to provide a quick primer on how <a>Kratos</a> code should be documented in-source, using doxygen.

Doxygen is the most popular tool for generating documentation from annotated source files of major languages (such as C++ or Python). It performs static code analysis on the repository to deduce a lot useful information from namespace contents to class composition and inheritance hierarchies. Although Doxygen can do this automatically without any comments or annotations, the proper description of classes and functions is essential for a useful documentation, which in turn is the cornerstone of attracting new users and developers to the project (and remembering what you were thinking a month ago when you wrote that class template).

## Documenting C++ Code

Doxygen is compatible with a vast variety of languages and the basic structure of annotating the code does not change, though minor syntax variations can occur between them. This section introduces essential tags and commands used throughout <a>Kratos</a>, with examples in C++.

Doxygen annotations are distinguished from regular comments by appending an extra comment character (```/``` or ```*```) to the comment signal (```//``` or ```/*``` respectively). Consecutive lines with annotations are considered to be part of the same annotation block, but regular comments without the exta characters are not included in the generated documentation.
```cpp
// A regular comment ignored by doxygen.

/// Single-line doxygen annotation.

/** Scoped doxygen annotation */

/// Another single-line annotation.
/** Followed by a multiline
 *  scoped doxygen annotation,
 *  in the same block.
 */
```

In C++, annotations should **directly precede** the constructs they document (**variable declarations**, **function declarations** or **class definitions**), and as such, should mostly be in header files. In case of forward-declared classes, annotations should only be included once, in front of the class definition.
```cpp
// A forward declaration of Node. Don't annotate this.
struct Node;

/** @brief Print data stored in a @ref Node to an output stream.
 *
 *  @param rStream: output stream to print the data to.
 *  @param rNode: node to extract data from.
 */
void PrintNode(std::ostream& rStream, const Node& rNode);

// Definition of Node => annotate it!
/// Simple 2D node, storing only its position and ID.
struct Node
{
    ///@name Data members
    ///@{

    /// Position on the x-axis.
    double m_x;

    /// Position on the y-axis.
    double m_y;

    ///@}

    ///@name Miscellaneous members
    ///@{

    /// Unique identifier of the node.
    unsigned int m_id;

    ///@}
};
```

[Here](data/example_1/html/index.html)'s what the generated documentation looks like for the example above. Notice the tags identified with <i>@</i> characters in the annotations; these are doxygen commands that relay formatting directives. They always modify the appearance of the text that follows them, but the number of characters/words they effect varies between commands. For example, the ```@ref``` command creates a link to the documentation of the name that follows it (```@ref Node``` links to ```Node``` in the generated documentation), while ```@param``` describes an argument of the function, taking a parameter name and a possibly multi-line description. There's also ```@name``` that tags every construct in its scope enclosed between ```@{``` and ```@}```. You can browse through the [complete list of doxygen commands](https://www.doxygen.nl/manual/commands.html) to get an idea of what you can work with, but here's a short list often used in <i>Kratos</i>:
- ```@brief```: brief description of the class/function/variable/concept that always gets displayed next to its name.
- ```@details```: detailed description of the class/function/variable/concept that is only displayed on the construct's own documentation page.
- ```@addtogroup```: add all constructs defined in the scope to a doxygen *module* (each *Kratos* application has its own *module*).
- ```@name```: Tag all constructs defined in the scope with a name, creating a separate paragraph for the in the documentation page. This command is used for member constructs within class definitions.
- ```@f(``` and ```@f)```: in-line latex equation.
- ```@f[``` and ```@f]```: latex equation on a new line.
- ```@param```: function argument description.
- ```@tparam```: template argument description
- ```@return```: return value description
- ```@ref```: create a link in the documentation to the name that follows it.
- ```@p```: reference to an argument local to the current annotation block.
- ```@a```: render the following name in italic.
- ```@b```: render the following name in bold.
- ```@c```: render the following name in typewriter font.


Let's take a deep dive and look at a reasonably well-decorated piece of code with lots of doxygen commands. Be sure to compare the source and the [generated documentation](data/example_2/html/index.html)
```cpp
/** @defgroup CompileTimeApplication
 *  @{
 *      @page Kratos Compile-time Application
 *      Utilities for performing calculations and logic
 *      at compile time to help with template programming.
 *  @}
 */

namespace Kratos {

///@addtogroup CompileTimeApplication
///@{

/** @brief Raise any base to an integer power.
 *
 *  @details Compute \f( a^n \f).
 *
 *  @tparam TNumeric: base type; can be any integer or floating point type (example: @c int or @c double).
 *  @param base: base to be raised to a power.
 *  @param exponent: multiply @p base by itself this many times
 *  @note this function can be invoked at compile time,
 *        and is guaranteed no to throw an exception.
 */
template <typename TNumeric>
constexpr TNumeric IntegerPower(TNumeric base, std::size_t exponent) noexcept
{
    if (exponent == 0) {
        return 1;
    }
    else {
        auto number_of_mults = 1;
        TNumeric power = base;
        while (number_of_mults < exponent) {
            if ((exponent - number_of_mults) / number_of_mults) {
                // Enough room to square the current state.
                power *= power;
                number_of_mults *= 2;
            }
            else {
                // Cannot square the current state
                // => multiply once by base
                power *= base;
                ++number_of_mults;
            }
        } // while

        return power;
    } // exponent != 0
} // IntegerPower()


/// @brief Compute the sum of the first @a n terms of a geometric series.
///
/// @details Compute the following expression: \f[ c \cdot \sum_{k=0}^n q^k \f]
///
/// @param coefficient: (@a c) coefficient multiplying each term of the series.
/// @param base: (@a q) base of the geometric series.
/// @param max_terms: (@a n) number of terms to sum up.
/// @note the classic reduced formula is used to compute the sum, unless @c base is @c 1.
///       Power calculation is deferred to @ref IntegerPower<double> which is less efficient
///       than @c std::pow but can be invoked at compile time and does not throw exceptions.
/// @warning be wary of floating point overflows.
constexpr double ComputeGeometricSum(double coefficient, double base, std::size_t max_terms) noexcept
{
    return coefficient * (base == 1 ? base * max_terms : (1 - IntegerPower(base, max_terms)) / (1 - base));
}

///@}

} // namespace Kratos
```

### Advice:
- Don't put annotations **inside function definitions**. Document everything in one block preceding the declaration.
- Define ```@addtogroup``` **inside namespaces** not outside of them, otherwise the constructs will be added to the documentation of the namespace but not the group.
- Always provide a brief description, and tag it with ```@brief```. It's very useful in the documentation.
- When using ```@ref``` to link to constructs in another namespace, make sure to specify the namespace as well. Doxygen won't know which construct to link to if there are more of them with the same name in different namespaces (or different languages - think of C++ classes and their bindings in python).
- If you're writing functions/classes that should not be in the documentation (helpers, implementation details), enclose them in ```namespace``` ```detail``` or ```impl``` if they're in headers, or [unnamed namespaces](https://en.cppreference.com/w/cpp/language/namespace#Unnamed_namespaces) if they're in source files.


## Documenting Python Code

For the most part, annotating python code is identical to documenting C++, save for the different comment character ```#``` and different languages features (such as package scopes instead of ```namespace```s). However, python is dynamically typed so type information must be manually provided if possible. Thankfully, python already provides a built-in remedy that is useful during development as well: [type hinting](https://peps.python.org/pep-0484/).

### Type Hints
Kratos developers seem to be oblivious to this feature of the language, so here's quick intro:
```py
# Type hinting a variable:
an_integer: int = 0 # Tell the interpreter that 'an_integer' is supposed to be an int

# Type hinting a function:
# The arguments' types can be hinted just like regular variables'
# The name followed by the arrow '->' after the argument list denotes the return type of the function.
def GetTime(model_part: KratosMultiphysics.ModelPart) -> float:
    return model_part.ProcessInfo(KratosMultiphysics.TIME)

# Type hinting optional arguments in a function:
def NewModelPart(model: KratosMultiphysics.Model, name: str = "Main") -> KratosMultiphysics.ModelPart:
    return model.CreateModelPart(name)

# Compound types (available since python3.9):
list_of_processes: list[KratosMultiphysics.Process] = [] # the interpreter will assume this list holds processes
mixed_tuple: tuple[str, int, float] = ("some_string", 1, 1.0)
material_parameters: dict[str, float] = {"density" : 2700.0, "youngs_modulus" : 6.9e10}

# Use the built-in 'typing' module for more complex hinting (available since python3.8)
import typing
cosine: typing.Callable[[float], float] = lambda x: math.cos(x)

# If the hinted type is not available in the script, or you are using an old version
# of python that doesn't support array hints, you can hint the types in strings instead.
# However, this won't be too useful for your IDE or the documentation, so this should be
# a last resort.
def GetOrigin() -> "tuple[float]":
    return (0.0, 0.0, 0.0)
```

Type hinted python code has the added benefit that your IDE will be able to provide you with more support (such as auto-complete suggestions).

### Annotations

As for annotations, the main difference is that functions' and classes' documentation must be written **inside** their **docstrings**. Variables' annotations are the same as in C++ and must immediately precede their definition. Python has no exact equivalent feature to C++'s namespace, so excluding parts of the code must be done manually. The best option to do this is to enclose the undocumented region between ```@cond impl``` and ```@endcond``` tags (the name of the conditional region in this case is ```impl``` but you are free to choose something else).

*The docstrings must not contain blank lines, otherwise doxygen commands won't be parsed in it. This is probably a doxygen bug or a configuration error on my part (let me know if you find out which).*

Here's an example and its [generated documentation](data/example_3/html/index.html) to give you an idea what python annotations should look like in practice:
```py
import typing

##@addtogroup utilities
##@{

class Memoizer:
    """@brief Utility class for caching and recalling results of expensive calculations.
    @details The memoized function must be an instance of a class implementing the @a __call__
    method. Unfortunately, normal functions are immutable in Python and so cannot be replaced.
    """

    def __init__(self, function: typing.Callable[[int], int]):
        """@brief Construct a memoizer from a callable object.
        @param function: Callable object to be memoized.
        @warning @p function must not keep internal state. Otherwise it won't
        always reproduce the same result and memoizing it won't make any sense.
        """
        self.__function: typing.Callable[[int], int] = function
        self.__cache: dict[int, int] = {}
        self.__backup: typing.Callable[[int], int] = function.__call__

    def __enter__(self) -> "Memoizer":
        """@brief Begin rerouting invocations."""
        type(self.__function).__call__ = lambda instance, argument: Memoizer.__call__(self, argument)
        return self

    def __call__(self, argument: int) -> int:
        """@brief Perform a lookup in the internal cache and invoke the wrapped function if that fails.
        @details Calls to the original function are redirected here first.
        @param argument: argument of the wrapped @a function.
        @return identical to the return value of the wrapped @a function
        """
        # Check whether this call is cached
        value: typing.Optional[int] = self.__cache.get(argument, None)
        if value == None: # It isn't => dispatch the original function and cache the result!
            value = self.__backup(argument)
            self.__cache[argument] = value
        return value

    def __exit__(self, *args) -> None:
        """@brief Stop rerouting invocations."""
        # For the sake of brevity, no error handling is done in this example.
        self.__function.__call__ = self.__backup


def Memoize(function: typing.Callable[[int], int]) -> Memoizer:
    """@brief Construct a @ref Memoizer for handling memoized contexts.
    @param function: callable to be memoized in the scope of the context.
    """
    return Memoizer(function)
##@}


##@cond impl
class FibonacciImpl:
    def __call__(self, n: int) -> int:
        return 0 if n < 1 else 1 if n==1 else self(n - 2) + self(n - 1)
##@endcond


##@defgroup maths

def Fibonacci(n: int, memoizer: Memoizer = Memoizer(FibonacciImpl())) -> int:
    """@brief Compute the @a n-th term of the fibonacci series.
    @details @f[
        F(n) = F(n-1) + F(n-2) \\
        F(0) = 0 \\
        F(1) = 1
    @f]
    @param n: 0-based index of the term to be computed.
    @param memoizer: instance of a @ref Memoizer for avoiding runaway recursion.
    @warning Term indices begin with 0!
    @ingroup maths
    """
    # A class instance is needed to overwrite __call__
    fibonacci = FibonacciImpl()
    with memoizer as MemoizedFibonacci:
        return MemoizedFibonacci(n)
```
