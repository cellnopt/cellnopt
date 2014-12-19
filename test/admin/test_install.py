from cno.admin import install



def test_install():
    install.install_all_cellnopt_dependencies(verbose=True)
    install.__main__()
