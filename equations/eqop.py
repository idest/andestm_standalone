#def applyfunc(self, op):
#    """Apply the same function to both sides of the equation.
#    Calls the function for each side, passing the lhs/rhs as
#    an argument, and then returns a new equation consisting of
#    the returned values."""
#    return self.func(op(self.lhs), op(self.rhs))
#
#def applylhs(self, op):
#    """Call the function on the lhs and return a new equation
#    whose lhs is the result of that function and rhs is the same
#    as in the original equation."""
#    return self.func(op(self.lhs), self.rhs)
#
#def applyrhs(self, op):
#    """Call the function on the rhs and return a new equation
#    whose rhs is the result of that function and lhs is the same
#    as in the original equation."""
#    return self.func(self.lhs, op(self.rhs))

def bsop(self, op):
    """Abbreviation for Both Sides Operation.
    Apply the same function to both sides of the equation.
    Calls the function for each side, passing the lhs/rhs as
    an argument, and then returns a new equation consisting of
    the returned values."""
    return self.func(op(self.lhs), op(self.rhs))

def lhsop(self, op):
    """Abbreviation for Left Hand Side Operation.
    Call the function on the lhs and return a new equation
    whose lhs is the result of that function and rhs is the same
    as in the original equation."""
    return self.func(op(self.lhs), self.rhs)

def rhsop(self, op):
    """Abbreviation for Right Hand Side Operation.
    Call the function on the rhs and return a new equation
    whose rhs is the result of that function and lhs is the same
    as in the original equation."""
    return self.func(self.lhs, op(self.rhs))
from sympy.core.relational import Relational
Relational.bsop = bsop
Relational.lhsop = lhsop
Relational.rhsop = rhsop

#test
#from sympy import Eq
#from sympy.abc import x, y
#eq = Eq(x,y)
#
#def op(x):
#    return x**2
#
#assert Eq(x, y).bsop(op) == Eq(x**2, y**2)
#assert Eq(x,y).bsop(lambda x: x*-1)) == Eq(-x, -y)
