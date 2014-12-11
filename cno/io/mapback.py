import networkx as nx
from cno.io import *
import pylab
from cno import CNOError



class MapBack(object):
    """

    :param pknmodel: compulsary and must be a CNOGraph data structure.
    :param model: a logical model, which shoudl be a subgraph of pknmodel with possible
        AND gates (see CNOGraph)
    :param links2map: the links present in model, which should be map back onto the pknmodel
        Could be a list reactions, a cnograph, or a Models (with one row).



    """

    def __init__(self, pknmodel, model):
        self.pknmodel = pknmodel
        self.model = model

    def search_path(self, source, target):
        """Search all simple paths in the graph G from source to target

        :param G:
        :param source:
        :param target:
        :return: list of edges
        """
        return list(nx.all_simple_paths(self.pknmodel, source, target))

    def plot_pknmodel(self, layout='inout'):
        self.pknmodel.plot(rank_method=layout)

    def plot_model(self, layout='inout'):
        self.model.plot(rank_method=layout)

    def plot_optimised_model(self, reactions):
        model = self.model.copy()
        mapback = [(reac,1) if reac in reactions else (reac,0)
                for reac in model.reactions]
        mapback = dict(mapback)
        model.set_edge_attribute('mapback', mapback)
        model.plot(edge_attribute='mapback', cmap='gray_r')
        return model

    def plot_mapback(self, reactions):
        """reactions should be the output of :meth:`mapback`.

        """
        model = self.pknmodel.copy()
        mapback = [(reac,1) if reac in reactions else (reac,0)
                for reac in model.reactions]
        mapback = dict(mapback)
        model.set_edge_attribute('mapback', mapback)
        model.plot(edge_attribute='mapback', cmap='gray_r')
        return model

    def plotall(self, reactions):
        pylab.figure(1)
        self.plot_pknmodel()
        pylab.title("pknmodel")

        pylab.figure(2)
        self.plot_model()
        pylab.title("model")

        pylab.figure(3)
        self.plot_optimised_model(reactions)
        pylab.title("optimised")

        reactions2 = self.mapback(reactions)
        pylab.figure(4)
        self.plot_mapback(reactions2)
        pylab.title("mapback")

    def mapback(self, links2map):
        """
        :param list links2map: list of edges 2 map from model to pknmodel.
            list of reactions

        :return: list of reactions in pknmodel that are on

        #the mapback for each link A->B in the compressed model is done looking
        #at the PKN as a graph and considering a subgraph of it including only
        #node A, node B and compressed nodes. All paths going from A to B in this
        #subnetwork are marked as active in the PKN.

        ::

            >>> c1 = CNOGraph(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyPB.csv"))
            >>> c2 = CNOGraph(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyPB.csv"))
            >>> c2.preprocessing()
            >>> newedges = mapback.mapback(c1,c2,edges)
            >>> c3 = CNOGraph(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyPB.csv"))
            >>> for edge in c3.edges():
            ...     if edge not in newedges:
            ...         c3.remove_edge(edge[0], edge[1])
            >>> c3.plotdot()

        .. todo:: current issue is nonc not being linked

        :References: although the implementation is slgihtly different, most of the code
            and testing were inspired by the mapBack function from CellNOptR 1.7 R package
            developed by F. Eduati.
        """
        for reaction in links2map:
            if reaction not in self.model.reactions:
                raise CNOError("Found an invalid reaction {0}".format(reaction) +
                        ", which is not in the model")

        compressed = [node for node in self.pknmodel.nodes()
            if node not in self.model.nodes()]

        # AND gates are handled as follows.
        # Case1, the AND gate is in the PKN. Nothing to do. The following algorithm shoud work
        # Case2, the AND gate is not in the PKN. In such case, we can just remove it
        #        and add the corresponding OR gates
        #
        #        If the AND gate is not in the PKN, it could be that it was renamed
        ands = [Reaction(x) for x in links2map if "^" in x]
        links2map = [x for x in links2map if "^" not in x]

        andsPKN = []
        for reaction in ands:
            if reaction.name in self.pknmodel.reactions:
                # no need to search for an AND gate that is in the PKN
                andsPKN.append(reaction.name)
            else:
                simple_reactions = reaction.ands2ors(reaction)
                links2map.extend(simple_reactions)

        newedges = []
        # for each link 2 map in the submodel
        additional_ands = []
        for reaction in links2map:
            reaction = Reaction(reaction)
            rhs = reaction.rhs
            lhs = reaction.get_signed_lhs_species()
            species = lhs['+'] + lhs['-']
            assert len(species) == len(list(set(species)))

            for startNode in species:
                # consider subgraph with startnode,end node and compressed nodes
                okNodes = compressed + [rhs, startNode]
                noNodes = [n for n in self.pknmodel.nodes() if n not in okNodes]
                #here is the subgrqph
                gg = self.pknmodel.copy()
                gg.remove_nodes_from(noNodes)
                #print(gg)

                # seqrch this sugrqph for paths between startnode and endnode
                paths = self.search_path(startNode, rhs)

                for path in paths:
                    for j2 in range(1, len((path))):
                        n1 = path[j2-1]
                        n2 = path[j2]
                        # FIXME iNOT are lost here
                        if startNode in lhs['+']:
                            reacPKN = n1 + "=" + n2
                        else:
                            reacPKN = "!" + n1 + "=" + n2
                        newedges.append(reacPKN)
                        if "^" in n1:
                            additional_ands.append(n1)
                        if "^" in n2:
                            additional_ands.append(n2)


        newedges = sorted(list(set(newedges)))
        newedges = [x for x in newedges if x.count("=")==1]
        # some edges may have been added from path that are actually not present
        # in the optimised networks. We have to identify the common edges in
        # the pkn and processed models. Then, to remove them in not part of the optimised
        # model.
        e1 = set(self.pknmodel.reactions)
        e2 = set(self.model.reactions)
        common_reactions = list(e1.intersection(e2))
        toremove = [r for r in common_reactions if r not in links2map]

        newedges = [edge for edge in newedges if edge not in toremove]
        newedges.extend(andsPKN)
        newedges.extend(list(set(additional_ands)))
        return newedges






