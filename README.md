# About the C++ code structure

The general coding style is to :
* use as few methods as possible
* use as many normal ("static") functions in methods as possible

This is why methods are generally written as a single function call, with many normal functions alongside them

This allows for several things :
* normal functions don't need object instanciation to be tested
* normal functions rarely use more arguments than needed (hence easier to test), and depend less on state
* normal functions can be reused and composed without inheritance or composition

One flaw of this approach is that some functions may have a very long list of arguments (we could abstract
the arguments in a "parameters" object, but this would lead to some instanciation overhead and less clarity)

Yet, we still use some object-oriented programming, the reasons are the following :
* C++ is object-oriented
* classes/structs are a way to store different functions arguments, and give them an additional consistency : 
for instance, when instantiating a BSEuroImplicit, which wraps an implicit Black Scholes PDE european pricing, we want to make
sure that it solves a PDE which has the same parameters and boundary conditions as a European Option
* inheritance can be a friendlier replacement for sum types ("union" in C, or the pipe syntax in ML-style languages). We use this to classify both EuropeanCallOption and EuropeanPutOption as EuropeanOption for instance (however it adds some overhead and makes refactoring more tedious)
* potential users of the library could choose between calling our pricing functions or instanciating objects and calling methods on these objects
* abstract classes give a template for other possible implementations