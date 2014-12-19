from cno import CNORfuzzy, cnodata



def test_fuzzy():
    c = CNORfuzzy(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
    c.optimise(2, maxgens=5)
    c.plot_mses()
    c.plot_errors()
    c.plot_fitness()

    score = c.simulate([ 3, 3, 3, 0, 4, 5, 0, 4, 5, 5, 6, 5, 2, 7, 6, 5, 0, 1, 0])
    assert score == 0.0254795668606642

    score = c.simulate(c.create_random_parameters())

    c.create_report()



