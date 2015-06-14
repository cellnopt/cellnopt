from biokit.rtools import RPackageManager


__all__ = ["install_all_cellnopt_dependencies"]


def install_all_cellnopt_dependencies(verbose=True):
    """script to install all dependencies for CellNOptR packages

        >>> from cno.admin import install_all_cellnopt_dependencies
        >>> install_all_cellnopt_dependencies()

    """
    pm = RPackageManager()
    installed_packages = pm.installed.index

    # requires those packages to be installed !!
    #packages = ['CellNOptR', 'CNORdt', 'CNORode', 'CNORfeeder', 'CNORfuzz']


    cellnopt = ["hash", "Rgraphviz", "RBGL", "graph", "RUnit", "igraph", "XML", "ggplot2", "RCurl"]
    ode= ["Rsolnp", "snowfall", "genalg"]
    feeder = ["catnet","minet"]
    meigor = ["Rsge"]
    fuzzy = ["xtable", "nloptr"]
    dt = []

    packages = sorted(list(set(cellnopt + meigor + fuzzy + feeder + ode + dt)))

    for package in packages:
        if package not in installed_packages:
            if verbose:
                print("Installing %s " % package)
            pm.install(package)
        else:
            if verbose:
                print("%s already installed. skipped" % package)

    # Rsge not maintained anymore so need to get it from arhive
    if "Rsge" not in installed_packages:
        pm.install_packages("http://cran.r-project.org/src/contrib/Archive/Rsge/Rsge_0.6.3.tar.gz")
    else:
        if verbose:
            print("%s already installed. skipped" % "Rsge")

    #MEIGOR
    pm.install("snowfall")
    pm.install("Rsolnp")
    if pm.is_installed("MEIGOR") is False:
        pm.install_packages("http://www.cellnopt.org/downloads/MEIGOR_0.99.6_svn3222.tar.gz",
            type="source")


if __name__ == "__main__":
    install_all_cellnopt_dependencies()
