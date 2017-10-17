############## Dedekind Domains Category #####

from sage.categories.integral_domains import IntegralDomains

class DedekindDomains(Category):
    
    def super_categories(self):
        return [IntegralDomains()]

    class ParentMethods:

        def ideal_with_rep(self, *gens, **kwds):
            '''
            Temporary ideal creation method. In general, the ideal constructor in sage is a bit
            complicated and warrents further attention.

            WARNING: This method does not check if the generators are recieved in the correct order.
            '''

            if len(gens) == 0:
                raise ValueError("needs at least one argument")
            if len(gens) > 2:
                raise ValueError("too many generators")
            if len(gens) == 1:
                return Ideal_dedekind(self, gens[0], gens[0], kwds)
            else:
                return Ideal_dedekind(self, gens[0], gens[1], kwds)

        def degree(self):
            '''
            This function should absolutely be removed, but is present for instructive purposes.
            '''
            return "override_after_category"

        def other_degree(self):
            return "the other degree"

    class ElementMethods:
        pass


def test_DedekindDomains_category():
    R.<t> = PolynomialRing(ZZ)
    OL = R.quo(t^2+1)

    # First lets check how things work by default
    print "Category: ", OL.category()
    print "degree output: ", OL.degree()  # We use this to check how override works
    
    try:
        print "New method test:", OL.ideal_with_rep()
    except AttributeError:
        print "OL does not have ideal_with_rep method"

    # Change the category. This tends to be done on an object initialization but constructor ``quo``
    # doesn't allow you to specify a category option

    OL._refine_category_(DedekindDomains())
    
    # Check how things work now
    print "Category: ", OL.category()
    print "degree output: ", OL.degree()  # Notice this is still the same!
    
    try:
        print "New method test: ", OL.other_degree() # now this works
    except AttributeError:
        print "OL does not have ideal_with_rep method"

    # if we really want to, we can actually assign methods
    OL.degree = OL.other_degree
    print "New degree output:", OL.degree()

    # TODO: declare various required methods needed to execute Dedekind domain algorithms.

class FieldsWithDedekindRingOfIntegers(Category):
    def super_categories(self):
        return [QuotientFields()]

