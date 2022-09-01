import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from gekko import GEKKO
from scipy.sparse import coo_matrix, eye, hstack


class Network:

    def __init__(self, path, disp=False):
        """
        In the initialization, the path of the xlsx file that contains the network
        data is loaded, the information is extracted from this same file, the
        incidence matrix is built, which describes the topology of the system and
        finally the upper and lower limits are established in each element of the
        system.
        """
        self.path = path
        self.disp = disp
        path = self.path

        # The information is extracted from the xlsx file
        self.node_info = pd.read_excel(path, sheet_name='node.info')
        self.node_info['node_id'] = self.node_info['node_id'].astype('int')
        self.node_info['type'] = self.node_info['type'].astype('int')

        self.node_dem = pd.read_excel(path, sheet_name='node.dem')
        self.node_dem['Total'] = self.node_dem.sum(axis=1)
        self.node_user = self.node_dem[self.node_dem['Total'] != 0]

        self.node_demcost = pd.read_excel(path, sheet_name='node.demcost')

        self.well = pd.read_excel(path, sheet_name='well')
        self.well['node'] = self.well['node'].astype('int')

        self.pipe = pd.read_excel(path, sheet_name='pipe')
        self.pipe['fnode'] = self.pipe['fnode'].astype('int')
        self.pipe['tnode'] = self.pipe['tnode'].astype('int')

        self.comp = pd.read_excel(path, sheet_name='comp')
        self.comp['fnode'] = self.comp['fnode'].astype('int')
        self.comp['tnode'] = self.comp['tnode'].astype('int')
        self.comp['Type'] = self.comp['Type'].astype('int')

        self.sto = pd.read_excel(path, sheet_name='sto')
        self.sto['node'] = self.sto['node'].astype('int')

        self.coordinates = pd.read_excel(path, sheet_name='coordinates')

        self.max_ratio = self.comp['ratio']
        self.N = len(self.node_info)
        self.W = len(self.well)
        self.P = len(self.pipe)
        self.C = len(self.comp)
        self.S = len(self.sto)
        self.Kij = self.pipe['Kij'].values
        N = len(self.node_info)
        # Incidence Matrix is created
        # Wells ok
        df_wells = pd.read_excel(path, sheet_name='well')
        W = len(self.well)
        wells = coo_matrix((np.ones(W, ), (self.well['node'] - 1, np.arange(W))), shape=(N, W))

        # Pipes ok
        df_pipes = pd.read_excel(path, sheet_name='pipe')

        P = len(self.pipe)
        data = np.concatenate((-1.0 * np.ones(P), np.ones(P)))
        row = pd.concat((self.pipe['fnode'] - 1, self.pipe['tnode'] - 1))
        col = np.concatenate(2 * [np.arange(P)])
        pipes = hstack(2 * [coo_matrix((data, (row, col)), shape=(N, P))]).toarray()
        # print('Pipes:', pipes.shape)

        # Compressors ok
        C = len(self.comp)
        data = np.concatenate((-1.0 * np.ones(C), np.ones(C)))
        row = pd.concat((self.comp['fnode'] - 1, self.comp['tnode'] - 1))
        col = np.concatenate(2 * [np.arange(C)])
        comps = coo_matrix((data, (row, col)), shape=(N, C))

        # Users ok
        users = hstack(len(self.node_demcost.T) * [eye(N)])

        # Storage
        self.sto = pd.read_excel(self.path, sheet_name='sto')
        S = len(self.sto)
        sto = coo_matrix((np.ones(S, ), (self.sto['node'] - 1, np.arange(S))), shape=(N, S))
        sto = hstack([sto, -1.0 * sto])

        self.Minc = (hstack((wells, pipes, comps, users, sto))).toarray()

        # The limits of the elementes of the network are stablished
        fp = self.pipe['FG_O'] * (self.pipe['FG_O'] > 0)
        fn = self.pipe['FG_O'] * (self.pipe['FG_O'] < 0)
        f_ext = np.concatenate((fp, fn))

        self.initial = np.concatenate((self.well['I'].values,
                                       f_ext,
                                       self.comp['fgc'].values,
                                       [0] * len(self.node_demcost) * len(self.node_demcost.T),
                                       [0] * len(self.sto.values),
                                       [0] * len(self.sto.values),
                                       self.node_info['p'].values))

        self.lb = np.concatenate((self.well['Imin'].values,
                                  [0] * len(self.pipe),
                                  self.pipe['Fg_min'].values,
                                  [0] * len(self.comp),
                                  [0] * len(self.node_dem) * (len(self.node_dem.T) - 1),
                                  [0] * len(self.sto) * 2,
                                  self.node_info['Pmin'].values))

        self.ub = np.concatenate((self.well['Imax'].values,
                                  self.pipe['Fg_max'].values,
                                  [0] * len(self.pipe),
                                  self.comp['fmaxc'],
                                  self.node_dem['Res'],
                                  self.node_dem['Ind'],
                                  self.node_dem['Com'],
                                  self.node_dem['NGV'],
                                  self.node_dem['Ref'],
                                  self.node_dem['Pet'],
                                  self.sto['V0'] - self.sto['V_max'].values,
                                  self.sto['Vmax'] - self.sto['V0'].values,
                                  self.node_info['Pmax'].values))

        self.objective_function()

    def objective_function(self):

        # Initialize Model
        self.m = GEKKO(remote=True)
        self.m.options.SOLVER = 3  # IPOPT Solver

        cost = np.concatenate((self.well['Cg'].values,
                               self.pipe['C_O'].values,
                               -1 * self.pipe['C_O'].values,
                               self.comp['costc'].values,
                               self.node_demcost['al_Res'].values,
                               self.node_demcost['al_Ind'].values,
                               self.node_demcost['al_Com'].values,
                               self.node_demcost['al_NGV'].values,
                               self.node_demcost['al_Ref'].values,
                               self.node_demcost['al_Pet'].values,
                               (self.sto['C_S+'] - self.sto['C_V']).values,
                               -1 * (self.sto['C_S-'] - self.sto['C_V']).values,
                               self.sto['C_V'].values))

        NV = (len(self.well) + 2 * len(self.pipe) + len(self.comp) +
              (len(self.node_demcost) * len(self.node_demcost.T)) +
              2 * len(self.sto) + len(self.node_info))

        self.NV_obj = NV - len(self.node_info)

        limits = (self.initial, self.lb, self.ub)

        # x = [self.m.Var(value=value,lb=lb,ub=ub) for value, lb, ub in zip(self.initial, self.lb, self.ub)]
        x = [self.m.Var(lb=lb, ub=ub) for lb, ub in zip(self.lb, self.ub)]
        self.X = np.array(x)

        for i, value in enumerate(np.array(self.initial)):
            self.m.fix_initial(x[i], value)

        self.X = np.array(x)

        B_lb = np.array([0] * 2 * len(self.pipe))
        B_ub = np.array([1] * 2 * len(self.pipe))

        B = [self.m.Var(lb=lb, ub=ub, integer=True) for lb, ub in zip(B_lb, B_ub)]
        self.B = np.array(B)

        # for x, x0 in zip(self.X, self.initial):
        #   self.m.fix_initial(x, x0)

        self.J = cost.T @ np.concatenate((self.X[:self.NV_obj],
                                          self.sto['V0'].values))

        self.constraints()
        self.solve_network()

    def constraints(self, ):

        self.gas_balance()
        self.comp_ratio()
        self.weymouth()

    # (2.24)
    def gas_balance(self):

        ftrans_max = self.pipe['Fg_max'].values
        pipe_indx1 = len(self.well)
        pipe_indx2 = len(self.well) + len(self.pipe)

        # self.m.Equations([self.X[pipe_indx1+i] + self.X[pipe_indx2+i] <= j for i, j in enumerate(ftrans_max)])
        self.m.Equations([self.net_flow(self.X)[i] <= j for i, j in enumerate(ftrans_max)])

        loads = self.node_dem['Total'].values
        self.A = self.X[:self.NV_obj] @ self.Minc.T

        for i in range(len(self.node_dem)):
            load = self.m.Param(value=loads[i])
            self.m.Equation(self.A[i] == load)

    def comp_ratio(self):
        press = self.X[-self.N:]
        subA = self.Minc[:, self.W + 2 * self.P:self.W + 2 * self.P + self.C]
        g = np.prod(press.reshape(-1, 1) ** subA, 0)
        g2 = self.m.Param(value=1)

        for i, constrain in enumerate(g):
            g1 = self.m.Param(value=self.max_ratio[i])
            self.m.Equation(constrain <= g1)
            self.m.Equation(constrain >= g2)

    def net_flow(self, x):
        flow = x[self.W:self.W + self.P] + x[self.W + self.P:self.W + 2 * self.P]
        return flow

    def weymouth(self):
        f = self.net_flow(self.X).reshape(-1, )
        press = self.X[-self.N:]
        Minc_P = self.Minc[:, self.W:self.W + self.P]
        Minc_P_i = Minc_P.copy()
        Minc_P_j = Minc_P.copy()
        Minc_P_i[Minc_P_i > 0] = 0
        Minc_P_i = -1 * Minc_P_i
        Minc_P_j[Minc_P_j < 0] = 0
        p_min = self.node_info['Pmin']
        p_max = self.node_info['Pmax']
        p_min_i = p_min @ (Minc_P_i)
        p_max_j = p_max @ (Minc_P_j)
        p_max_i = p_max @ (Minc_P_i)
        p_min_j = p_min @ (Minc_P_j)
        phi_plus_lowlimit = p_min @ np.abs(Minc_P)
        phi_plus_uplimit = p_max @ np.abs(Minc_P)
        phi_minus_lowlimit = p_min_i - p_max_j
        phi_minus_uplimit = p_max_i - p_min_j

        phi_plus = (press @ np.abs(Minc_P))
        phi_minus = (press @ (-1 * Minc_P))

        y = [self.m.Var(integer=True) for i in range(self.P)]
        Phi = [self.m.Var() for i in range(self.P)]
        M = [self.m.Var(lb=1) for i in range(self.P)]

        Phi_lowlimit1 = (phi_plus_lowlimit * phi_minus) + \
                        (phi_plus * phi_minus_lowlimit) - \
                        (phi_plus_lowlimit * phi_minus_lowlimit)

        Phi_lowlimit2 = (phi_plus_uplimit * phi_minus) + \
                        (phi_minus * phi_minus_uplimit) - \
                        (phi_plus_uplimit * phi_minus_uplimit)

        Phi_uplimit1 = (phi_plus_uplimit * phi_minus) + \
                       (phi_plus * phi_minus_lowlimit) - \
                       (phi_plus_uplimit * phi_minus_lowlimit)

        Phi_uplimit2 = (phi_plus_lowlimit * phi_minus) + \
                       (phi_plus * phi_minus_uplimit) - \
                       (phi_plus_lowlimit * phi_minus_uplimit)

        for i in range(self.P):
            self.m.Equation(phi_plus[i] <= self.m.Param(value=phi_plus_uplimit[i]))
            self.m.Equation(phi_plus[i] >= self.m.Param(value=phi_plus_lowlimit[i]))
            self.m.Equation(phi_minus[i] <= self.m.Param(value=phi_minus_uplimit[i]))
            self.m.Equation(phi_minus[i] >= self.m.Param(value=phi_minus_lowlimit[i]))

            self.m.Equation(Phi[i] <= self.m.Param(value=Phi_uplimit1[i]))
            self.m.Equation(Phi[i] <= self.m.Param(value=Phi_uplimit2[i]))
            self.m.Equation(Phi[i] >= self.m.Param(value=Phi_lowlimit1[i]))
            self.m.Equation(Phi[i] >= self.m.Param(value=Phi_lowlimit2[i]))

            self.m.Equation(self.m.Param(value=f[i]) <= self.Kij[i] * Phi[i] + \
                            M[i] ** 2 * (1 - y[i]))

            self.m.Equation(self.m.Param(value=f[i]) <= -1 * self.Kij[i] * Phi[i] + \
                            M[i] ** 2 * y[i])

    def solve_network(self):
        self.m.options.OTOL = 1e-7
        self.m.options.MAX_ITER = 1e9
        self.m.solver_options = ['mu_strategy adaptive',
                                 'constr_viol_tol 1e-7',
                                 'acceptable_tol 1e-7', ]
        #  'bound_push 1e-10',\
        #  'bound_frac 1e-10']
        self.m.Minimize(self.J)
        self.m.open_folder()
        self.m.solve(disp=self.disp)

    def show_values(self):

        print('Objective: ' + str(self.m.options.objfcnval))
        for i in range(len(self.well)):
            inj = self.well['node'].values[i]
            print('Inj', inj, ':', self.X[i].value[0])

        print('\n')
        f_plus = self.X[self.W: self.W + self.P]
        f_minus = self.X[self.W + self.P: self.W + self.P + self.P]
        for i in range(self.P):
            From = self.pipe['fnode'][i]
            To = self.pipe['tnode'][i]
            print('Net flow', From, To, ':', round(f_plus[i][0] + f_minus[i][0], 3))

        print('\n')
        for i in range(len(self.comp)):
            From = self.comp['fnode'][i]
            To = self.comp['tnode'][i]
            print('comp', From, To, ':', self.X[len(self.well) + 2 * len(self.pipe) + i].value[0])

        print('\n')
        print('Objective: ' + str(self.m.options.objfcnval))

    def show_network(self, flow_max=None):
        nodes_ = np.arange(1, len(self.node_info) + 1)
        edges_pipe = self.pipe[['fnode', 'tnode']].values
        edges_comp = self.comp[['fnode', 'tnode']].values

        edges = np.concatenate((edges_pipe,
                                edges_comp))

        f_plus = self.X[self.W: self.W + self.P]
        f_minus = self.X[self.W + self.P: self.W + self.P + self.P]
        f_pipe = np.round([f_plus[i][0] + f_minus[i][0] for i in range(self.P)], 3)
        f_comp_ = self.X[self.W + 2 * self.P:self.W + 2 * self.P + self.C]
        f_comp = np.round([i[0] for i in f_comp_], 3)

        f = np.concatenate((f_pipe, f_comp))
        data = {'from': edges[:, 0],
                'to': edges[:, 1],
                'flow': f}
        graph = pd.DataFrame(data).sort_values('from')
        cmap = plt.cm.coolwarm
        fmax = f.max()
        if not flow_max:
            fmax = f.max()

        elif (type(flow_max) is float or
              type(flow_max) is int):
            fmax = np.abs(flow_max)
        fmin = -fmax
        # G = nx.from_pandas_edgelist(graph, 'from', 'to', create_using=nx.Graph() )
        nodes = self.coordinates
        pos = {nodes['id'][i]: [nodes['x'][i], nodes['y'][i]] for i in range(len(nodes))}
        G = nx.DiGraph()
        G.add_nodes_from(nodes_)
        G.add_edges_from(edges)
        options = {'node_size': 300,
                   'font_size': 9,
                   'width': 3,
                   'node_shape': '.',
                   'with_labels': True,
                   'node_color': 'orange',
                   'font_weight': 'normal',
                   'edge_color': graph['flow'],
                   'edge_vmin': fmin,
                   'edge_vmax': fmax,
                   'edge_cmap': cmap,
                   'pos': pos
                   }

        fig = plt.figure(dpi=250)
        # nx.draw_networkx(G, )
        nx.draw(G, **options)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=fmin, vmax=fmax))
        sm.set_array([])
        cbar = plt.colorbar(sm)
        plt.show()


if __name__ == "__main__":
    path = '/home/cristian/workspace/PycharmProjects/natural_gas/estructura_red_2022_base.xlsx'
    p = Network(path, disp=True)
    p.show_values()
    p.show_network(flow_max=None)
