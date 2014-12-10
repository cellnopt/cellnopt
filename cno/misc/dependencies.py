import easydev
from cno import CNOGraph


def plot_dependencies(package='cno', show=False, filename=None):

    main = easydev.dependencies.get_dependencies(package)

    # first, we fetch all dependencies
    c = CNOGraph()
    deps = easydev.dependencies.get_dependencies(package)
    package_version = [dep.version for dep in deps if dep.project_name ==  package][0]
    for dep in deps:
        version = dep.version
        name= dep.project_name
        c.add_reaction(package+'-'+package_version + "="+name + "-" + version)

    #actually, dependencies themeselves depends on other packages
    # that are all present in deps variable but we are missing 
    # the overall DAG so let us loop over the dependencies
    deps = [dep for dep in deps if
            len(easydev.dependencies.get_dependencies(dep.project_name))>1]

    newdeps = []
    count = 0
    while len(deps) > 1 and count < 2:
        deps = [dep for dep in deps if
                len(easydev.dependencies.get_dependencies(dep.project_name))>1]
        for dep in deps:
            for this in easydev.dependencies.get_dependencies(dep.project_name):
                if this.project_name != dep.project_name:
                    newdeps.append(this)
                version = dep.version
                name= dep.project_name
                c.add_reaction(this.project_name +"-" +this.version+"="+name + "-" + version)
        deps = newdeps
        count +=1

    c.remove_self_loops()

    # keep only longest path. Not use this is the good way of doing it...
    c2 = CNOGraph()

    import networkx as nx
    for clique in list(nx.find_cliques(c.to_undirected())):
        for i in range(0, len(clique)-1):
            c2.add_edge(clique[i], clique[i+1])
    c2.plot(filename=filename, show=show)




