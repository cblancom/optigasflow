import pandas as pd
import numpy as np 
from scipy.sparse import coo_matrix, eye, hstack
from gekko import GEKKO
import matplotlib.pyplot as plt
import dill

class InterconnectedEnergySystemOptimizer:
    def __init__(self, data_gas, data_power, T=1, disp=True, seed=42, 
                 num_approximations=20, weymouth_type='MPCC'):
        
        self.data_gas = data_gas
        self.data_power = data_power
        self.T = T
        self.disp = disp 
        self.seed = seed
        np.random.seed(self.seed)
        self.num_approximations = num_approximations
        self.weymouth_type = weymouth_type
        self.flag = False
        
        self.extract_data_gas(self.data_gas)
        self.extract_data_power(self.data_power)
        self.m = GEKKO(remote=True)
        self.m.options.SOLVER = 3
        
        self.objective_function_gas()
        self.objective_function_power()
        
        self.interconnection()
        self.gas_constraints()
        self.power_constraints()
        self.solve_network()

    def interconnection(self):
        self.gas_to_power()
        # self.power_to_gas()

    def gas_constraints(self):
        self.gas_balance()
        self.comp_ratio()
        if self.weymouth_type == 'MPCC':
            self.weymouth_MPCC()
        elif self.weymouth_type == 'SOC':
            self.weymouth_SOC()
        elif self.weymouth_type == 'Taylor':
            self.weymouth_Taylor()
        self.storage_constraint()

    def power_constraints(self):
        self.power_balance()
        self.line_power_flow()

    def extract_data_gas(self, data):
        self.node_info = pd.read_excel(data, sheet_name='node.info')
        self.node_info['node_id'] = self.node_info['node_id'].astype('int')
        self.node_info['type'] = self.node_info['type'].astype('int')

        self.node_dem = pd.read_excel(data, sheet_name='node.dem')
        try:
            self.node_dem = self.node_dem.drop(columns='Nodo')
        except:
            pass
        self.node_dem['Total'] = self.node_dem.sum(axis=1)
        self.node_user = self.node_dem[self.node_dem['Total'] != 0]
        self.loads_gas = self.node_dem['Total'].values
        
        self.node_demcost = pd.read_excel(data, sheet_name='node.demcost')

        self.well = pd.read_excel(data, sheet_name='well')
        self.well['node'] = self.well['node'].astype('int')

        self.node_deminitial = pd.read_excel(data, sheet_name='node.deminitial')

        self.pipe = pd.read_excel(data, sheet_name='pipe')
        self.pipe['fnode'] = self.pipe['fnode'].astype('int')
        self.pipe['tnode'] = self.pipe['tnode'].astype('int')

        self.comp = pd.read_excel(data, sheet_name='comp')
        self.comp['fnode'] = self.comp['fnode'].astype('int')
        self.comp['tnode'] = self.comp['tnode'].astype('int')
        self.comp['Type'] = self.comp['Type'].astype('int')

        self.sto = pd.read_excel(data, sheet_name='sto')
        self.sto['node'] = self.sto['node'].astype('int')
        try:
            self.coordinates = pd.read_excel(data, sheet_name='coordinates')
        except:
            self.coordinates = False
        self.max_ratio = self.comp['ratio']
        self.N = len(self.node_info)
        self.W = len(self.well)
        self.P = len(self.pipe)
        self.C = len(self.comp)
        self.UNS = len(self.node_demcost) * len(self.node_demcost.T)
        self.S = len(self.sto)
        self.Kij = self.pipe['Kij'].values
        N = len(self.node_info)

        # df_wells = pd.read_excel(self.data, sheet_name='well')
        # W = len(self.well)
        wells = coo_matrix(
            (np.ones(self.W, ), (self.well['node'] - 1, np.arange(self.W))), shape=(self.N, self.W))
        # Pipes ok
        # df_pipes = pd.read_excel(self.data, sheet_name='pipe')

        # P = len(self.pipe)
        data = np.concatenate((-1.0 * np.ones(self.P), np.ones(self.P)))
        row = pd.concat((self.pipe['fnode'] - 1, self.pipe['tnode'] - 1))
        col = np.concatenate(2 * [np.arange(self.P)])

        pipes = hstack(
            2 * [coo_matrix((data, (row, col)), shape=(self.N, self.P))]).toarray()
        # print('Pipes:', pipes.shape)
        # Compressors ok
        # C = len(self.comp)
        data = np.concatenate((-1.0 * np.ones(self.C), np.ones(self.C)))
        row = pd.concat((self.comp['fnode'] - 1, self.comp['tnode'] - 1))
        col = np.concatenate(2 * [np.arange(self.C)])
        comps = coo_matrix((data, (row, col)), shape=(N, self.C))
        # Users ok
        users = hstack(len(self.node_demcost.T) * [eye(N)])
        # Storage
        # self.sto = pd.read_excel(self.data, sheet_name='sto')
        # S = len(self.sto)
        sto = coo_matrix(
            (np.ones(self.S, ), (self.sto['node'] - 1, np.arange(self.S))), shape=(self.N, self.S))
        sto = hstack([sto, -1.0 * sto])
        self.Minc = (hstack((wells, pipes, comps, users, sto))).toarray()
        
        # print(np.array(self.node_deminitial).flatten(order='F').shape)
        
        self.NV_gas = self.W + 2*self.P + self.C + self.UNS + 3*self.S + self.N
        self.NV_obj_gas = self.NV_gas - self.N
        

    def extract_data_power(self,data):
        self.bus = pd.read_excel(data, sheet_name='bus')
        self.gen = pd.read_excel(data, sheet_name='gen')
        self.gencost = pd.read_excel(data, sheet_name='gencost')
        self.demand = pd.read_excel(data, sheet_name='demand')
        self.branch = pd.read_excel(data, sheet_name='branch')
        self.G2P = pd.read_excel(data, sheet_name='G2P')
        
        if 'Pmax' not in self.branch.columns:
            self.branch['Pmax'] = self.branch['angmax'] * np.pi/(180*self.branch['x'])
            self.branch['Pmin'] = self.branch['angmin'] * np.pi/(180*self.branch['x'])
        # self.loads = self.demand['load']
        self.loads_power = np.array(self.demand)[:, 0:self.T]
        if np.array(self.loads_power).shape[1] > self.T:
            raise ValueError("Insufficient demands per day")
        self.x_param = self.branch['x']
        # Num elements
        self.Nbus = len(self.bus)
        self.Ngen = len(self.gen)
        self.Nbranch = len(self.branch)
        self.Nload = len(self.demand)
        # Incidence Matrix
        self.bus_inc = coo_matrix((np.ones(self.Ngen, ), 
                                   (self.gen['bus'] - 1, np.arange(self.Ngen ))), 
                                  shape=(self.Nbus , self.Ngen ))
        data = np.concatenate((-1.0 * np.ones(self.Nbranch), 
                               np.ones(self.Nbranch)))
        row = pd.concat((self.branch['fbus'] - 1, 
                         self.branch['tbus'] - 1))
        col = np.concatenate(2 * [np.arange(self.Nbranch)])
        self.branch_inc = coo_matrix((data, (row, col)), shape=(self.Nbus, self.Nbranch))
        self.loads_inc = eye(self.Nbus)
        self.ps_minc = (hstack((self.bus_inc, self.loads_inc, self.branch_inc))).toarray() 
        self.NV_ps = self.Ngen + self.Nload + self.Nbranch + self.Nbus
        self.NV_obj_ps = self.NV_ps - self.Nbranch - self.Nbus

    def objective_function_gas(self):
        X = self.m.Array(self.m.Var, (self.T, self.NV_gas))
        param_cost = np.concatenate((self.well['Cg'].values,
                                     self.pipe['C_O'].values,
                                     -1 * self.pipe['C_O'].values,
                                     self.comp['costc'].values,
                                     self.node_demcost['al_Res'].values,
                                     self.node_demcost['al_Ind'].values,
                                     self.node_demcost['al_Com'].values,
                                     self.node_demcost['al_NGV'].values,
                                     self.node_demcost['al_Ref'].values,
                                     self.node_demcost['al_Pet'].values,
                                     self.sto['C_S+'].values,
                                     self.sto['C_S-'].values,
                                     self.sto['C_V'].values))
        for i in range(self.T):
            # The limits of the elementes of the network are stablished
            fp = self.pipe['FG_O'] * (self.pipe['FG_O'] > 0)
            fn = self.pipe['FG_O'] * (self.pipe['FG_O'] < 0)
            f_ext = np.concatenate((fp, fn))
            initial = np.concatenate((self.well['I'].values,
                                        f_ext,
                                        self.comp['fgc'].values,
                                        np.array(self.node_deminitial).flatten(order='F'),
                                        self.sto['fs+'],
                                        self.sto['fs-'],
                                        self.sto['V0'].values,
                                        self.node_info['p'].values))
            lb = np.concatenate((self.well['Imin'].values,
                                    [0] * len(self.pipe),
                                    self.pipe['Fg_min'].values,
                                    [0] * len(self.comp),
                                    [0] * len(self.node_dem) *
                                    (len(self.node_dem.T) - 1),
                                    [0] * len(self.sto),
                                    [0] * len(self.sto),
                                    [0] * len(self.sto),
                                    self.node_info['Pmin'].values))
            ub = np.concatenate((self.well['Imax'].values,
                                    self.pipe['Fg_max'].values,
                                    [0] * len(self.pipe),
                                    self.comp['fmaxc'],
                                    self.node_dem['Res'],
                                    self.node_dem['Ind'],
                                    self.node_dem['Com'],
                                    self.node_dem['NGV'],
                                    self.node_dem['Ref'],
                                    self.node_dem['Pet'],
                                    self.sto['fsmax'].values,
                                    self.sto['fsmax'].values,
                                    self.sto['Vmax'].values,
                                    self.node_info['Pmax'].values))
            # print(initial.shape, lb.shape, ub.shape)
            if (len(ub) != len(lb) or (len(ub) != len(initial))):
                    raise ValueError("Problem with the boundaries.")
            self.limits_gas = (initial, lb, ub)
            for j, (value, lb, ub) in enumerate(zip(initial, lb, ub)):
                X[i, j].value = value
                X[i, j].lower = lb
                X[i, j].upper = ub    
        # np.random.seed(self.seed)
        # self.well_costs = np.random.normal(4500, 1500, self.T)
        self.cost_gas = np.tile(param_cost, (self.T))
        # for i in range(len(self.cost_gas)):
        #     if i % len(param_cost) == 0:
        #         well_index = i // len(param_cost)  # integer division to get index in well array
        #         self.cost_gas[i] = self.well_costs[well_index]
        self.X_gas = X[:, :self.NV_obj_gas]
        self.gas_flow = X[:, :self.NV_obj_gas-self.S] 
        self.press = X[:, -self.N:]
        self.storage = X[:, -self.N-3*self.S:-self.N]
        self.storage[-1, 1::3] = self.m.Param(value=0)
        self.storage[-1, 2::3] = self.m.Param(value=0)
        
        if (len(self.cost_gas) != len(self.X_gas.flatten())):
            raise ValueError("Problem with the objective function.")

        for i in range(len(self.cost_gas)):
            self.m.Obj(self.cost_gas[i] * self.X_gas.flatten()[i])

    def objective_function_power(self):
        X = self.m.Array(self.m.Var, (self.T, self.NV_ps))
        param_cost = np.concatenate((self.gencost['cost'],
                                     [1e6] * self.Nload
                                     ))
        for i in range(self.T):
            initial = np.concatenate((self.gen['Pg'],
                                       [0] * self.Nload,
                                       [0] * self.Nbranch,
                                       [0] * self.Nbus
                                       ))
            lb = np.concatenate((self.gen['Pmin'],
                                    [0] * self.Nbus,
                                    self.branch['Pmin'],  
                                    [-6.3] * self.Nbus    # 360 to rad                                
                                    ))
            ub = np.concatenate((self.gen['Pmax'],
                                    self.loads_power[:, i],
                                    self.branch['rateA'],
                                    [6.3] * self.Nbus
                                    ))
            if (len(ub) != len(lb) or (len(ub) != len(initial))):
                raise ValueError("Problem with the boundaries.")
            self.limits_gas = (initial, lb, ub)
            for j, (value, lb, ub) in enumerate(zip(initial, lb, ub)):
                X[i, j].value = value
                X[i, j].lower = lb
                X[i, j].upper = ub
        # self.cost = np.tile(self.param_cost, (self.T))
        self.cost_power = np.tile(param_cost, (self.T))
        self.X_power = X[:, :self.NV_obj_ps]
        self.power_flow = X[:, self.NV_obj_ps:self.NV_obj_ps+self.Nbranch]
        self.flow_minc = X[:, :self.NV_ps-self.Nbus]
        self.delta = X[:, -self.Nbus:]
        if (len(self.cost_power) != len(self.X_power.flatten())):
            raise ValueError("Problem with the objective function.")
        for i in range(len(self.cost_power)):
            self.m.Obj(self.cost_power[i] * self.X_power.flatten()[i])
    
    def gas_to_power(self):
        self.loads = self.m.Array(self.m.Var, (self.T, len(self.loads_gas)))
        for i in range(self.T):
            for j in range(len(self.loads_gas)):    
                self.loads[i, j].value = self.loads_gas[j]
                self.loads[i, j].lower = self.loads_gas[j]
                self.loads[i, j].upper = self.loads_gas[j]
        X_gen = self.X_power[:, :self.Ngen]
        idx = self.G2P['node'] - 1
        bus = self.G2P['bus'] - 1
        heat_rate = self.G2P['heat_rate']
        for i in range(len(idx)):
            self.loads[:, idx[i]] =  X_gen[:,bus[i]] * heat_rate[i]

    def power_to_gas(self):
        self.loads_power_tmp = self.m.Array(self.m.Var, (self.T, len(self.loads_power)))
        
        for i in range(self.T):
            for j in range(len(self.loads_power)):    
                self.loads_power_tmp[i, j].value = self.loads_power[j, i]
                self.loads_power_tmp[i, j].lower = self.loads_power[j, i]
                self.loads_power_tmp[i, j].upper = self.loads_power[j, i]

        self.phi1 = self.m.Var(value=45, lb=45)
        self.phi2 = self.m.Var(value=67.5, lb=67.5)
        self.loads_power_tmp[0][4] = self.phi1 
        self.loads_power_tmp[1][4] = self.phi2 
        
        flow_comp = self.gas_flow[0, self.W+2*self.P:self.W+2*self.P+self.C]
        # self.m.Equation(self.phi1 * self.press[0, 2]**(0.15) == 4.88 * flow_comp[1] * (self.press[0, 3]**(0.15) - self.press[0, 2]**(0.15)))
        self.m.Equation(self.phi1 * self.press[0, 2]**(0.236) == flow_comp[1] * (self.press[0, 3]**(0.236) - self.press[0, 2]**(0.236)))
        # self.m.Equation(self.phi1 == flow_comp[1] * self.m.sqrt((self.press[0, 3]**(0.236)) - 1))


    def gas_balance(self):
        # if not self.loads:
        if not hasattr(self, 'loads'):
            self.loads = self.m.Array(self.m.Var, (self.T, len(self.loads_gas)))
            for i in range(self.T):
                for j in range(len(self.loads_gas)):
                    self.loads[i, j].value = self.loads_gas[j]
                    self.loads[i, j].lower = self.loads_gas[j]
                    self.loads[i, j].upper = self.loads_gas[j]

        if self.T ==1:
            ftrans_max = self.pipe['Fg_max'].values
            self.m.Equations([self.net_flow(self.gas_flow[0])[i] <=
                            j for i, j in enumerate(ftrans_max)])
            # self.A = self.X[time, :self.NV_obj] @ self.Minc.T
            self.A_gas = self.gas_flow @ self.Minc.T
            # print(self.A.shape)
            for i in range(len(self.node_dem)):
                load = self.m.Param(value=self.loads_gas[i])
                self.m.Equation(self.A_gas[0][i] == self.loads[0, i])
                # self.m.Equation(self.A_gas[0][i] == load)
        elif self.T >= 1:
            ftrans_max = self.pipe['Fg_max'].values
            
            for time in range(self.T):
                self.m.Equations([self.net_flow(self.gas_flow[time])[i] <=
                                j for i, j in enumerate(ftrans_max)])
                # self.A = self.X[time, :self.NV_obj] @ self.Minc.T
                self.A_gas = self.gas_flow[time] @ self.Minc.T
                # print(self.A.shape)
                for i in range(len(self.node_dem)):
                    load = self.m.Param(value=self.loads_gas[i])
                    self.m.Equation(self.A_gas[i] == self.loads[time, i])
                        # if i == 6:
                        #     self.m.Equation(self.A_gas[i] == self.X_power[time, 2] * self.G2P['heat_rate'][0]) 
                        #     # self.m.Equation(self.A_gas[i] == self.X_power[time, 2] ) 
                        # else:        
                        #     load = self.m.Param(value=self.loads_gas[i])
                        #     self.m.Equation(self.A_gas[i] == self.loads[time, i])
                        #     # self.m.Equation(self.A_gas[i] == load)

    def power_balance(self):
        if self.T == 1:
            ptrans_max = self.branch['rateA'].values
            self.m.Equations([self.power_flow[0][i] <= j for i, j in enumerate(ptrans_max)])
            self.A_power = self.ps_minc @ self.flow_minc.T
            for i in range(self.Nload):
                load = self.m.Param(value=self.loads_power[i])
                self.m.Equation(self.A_power[i][0] == load)
        elif self.T >= 1:
            for time in range(self.T):
                ptrans_max = self.branch['rateA'].values
                self.m.Equations([self.power_flow[time][i] <= j for i, j in enumerate(ptrans_max)])
                # self.A = self.ps_minc @ self.flow_minc[time]
                self.A_power = self.flow_minc[time] @ self.ps_minc.T
                for i in range(self.Nload):
                    load = self.m.Param(value=self.loads_power[i][time])
                    # load = self.m.Param(value=self.loads_power_tmp[time][i])
                    self.m.Equation(self.A_power[i] == load)
    
    def comp_ratio(self):
        for time in range(self.T):
            # press = self.X[time, -self.N:]
            subA = self.Minc[:, self.W + 2 * self.P:self.W + 2 * self.P + self.C]
            g = np.prod(self.press[time].reshape(-1, 1) ** subA, 0)
            g2 = self.m.Param(value=1)
            for i, constrain in enumerate(g):
                g1 = self.m.Param(value=self.max_ratio[i])
                self.m.Equation(constrain <= g1)
                self.m.Equation(constrain >= g2)


    def weymouth_MPCC(self):
        print('MPCC')
        for time in range(self.T):
            f = (self.net_flow(self.gas_flow[time]).reshape(-1,))
            press2 = np.array(self.press[time]) ** 2
            for i in range(self.P):
                i_index = (self.pipe['fnode'] - 1)[i]
                j_index = (self.pipe['tnode'] - 1)[i]
                p = press2[i_index] - press2[j_index]
                K = self.Kij[i]
                self.m.Equation(self.m.sign2(f[i]) * f[i]**2 == p*K)

    def weymouth_SOC(self):
        print('SOC')
        for time in range(self.T):
            f = self.net_flow(self.gas_flow[time]).reshape(-1,)
            press = self.press[time]
            Minc_P = self.Minc[:, self.W:self.W+self.P]
            Minc_P_i = Minc_P.copy()
            Minc_P_j = Minc_P.copy()
            Minc_P_i[Minc_P_i > 0] = 0
            Minc_P_i = -1*Minc_P_i
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
            phi_minus = (press @ (-1*Minc_P))
            
            
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
                self.m.Equation(phi_plus[i] <= self.m.Param(value = phi_plus_uplimit[i]))
                self.m.Equation(phi_plus[i] >= self.m.Param(value = phi_plus_lowlimit[i]))
                self.m.Equation(phi_minus[i] <= self.m.Param(value = phi_minus_uplimit[i]))
                self.m.Equation(phi_minus[i] >= self.m.Param(value = phi_minus_lowlimit[i]))

                self.m.Equation(Phi[i] <= self.m.Param(value = Phi_uplimit1[i]))
                self.m.Equation(Phi[i] <= self.m.Param(value = Phi_uplimit2[i]))
                self.m.Equation(Phi[i] >= self.m.Param(value = Phi_lowlimit1[i]))
                self.m.Equation(Phi[i] >= self.m.Param(value = Phi_lowlimit2[i]))

                self.m.Equation(self.m.Param(value = f[i]) <= self.Kij[i] * Phi[i] + \
                                M[i]**2 * (1 - y[i]))
                
                self.m.Equation(self.m.Param(value = f[i]) <= -1*self.Kij[i] * Phi[i] + \
                                M[i]**2 * y[i])

    def weymouth_Taylor(self):
        print('Taylor')    
        for time in range(self.T):
            num_approximations = self.num_approximations
            press = self.press[time]
            f = self.net_flow(self.gas_flow[time]).reshape(-1,)

            for i in range(self.P):
                low_limit = self.node_info['Pmin'].values[i]
                upp_limit = self.node_info['Pmax'].values[i]
                K = self.Kij[i]
                pi = press[self.pipe['fnode'][i] - 1]
                pj = press[self.pipe['tnode'][i] - 1]
                g_value = self.m.Param(value = f[i])

            for j in range(num_approximations):
                PO, PI = np.sort(np.random.uniform(low=low_limit, high=upp_limit, size=(2,)))
                g = K*( (PI/np.sqrt(PI**2-PO**2))*pi + (PO/np.sqrt(PI**2-PO**2))*pj)
                self.m.Equation(self.m.abs(g_value) <= g)

    def line_power_flow(self,):
        if self.T == 1:
            f = (self.power_flow[0].reshape(-1,))
            delta = np.array(self.delta[0])
            for i in range(self.Nbranch):
                i_index = (self.branch['fbus'] - 1)[i]
                j_index = (self.branch['tbus'] - 1)[i]
                d = (delta[i_index] - delta[j_index]) / self.x_param[i]
                # K = self.Kij[i]
                self.m.Equation(f[i] == d)
        elif self.T >= 1:
            for time in range(self.T):
                f = (self.power_flow[time].reshape(-1,))
                delta = np.array(self.delta[time])
                for i in range(self.Nbranch):
                    i_index = (self.branch['fbus'] - 1)[i]
                    j_index = (self.branch['tbus'] - 1)[i]
                    d = (delta[i_index] - delta[j_index]) / self.x_param[i]
                    # K = self.Kij[i]
                    self.m.Equation(f[i] == d)

    def net_flow(self, x):
        flow = x[self.W:self.W + self.P] + \
            x[self.W + self.P:self.W + 2 * self.P]
        return flow

    def storage_constraint(self,):
            self.storage[0, -self.S:] = self.m.Param(value=self.sto['V0'].values)
            # self.storage[-1, :] = self.m.Param(value=0)
            f_plus = self.storage[:, :self.S]
            f_minus = self.storage[:, self.S:2*self.S]
            V = self.storage[:, -self.S:]
            # sto_var = self.X[:, -self.N-self.S*2:-self.N]
            if self.T == 1:
                for sto in range(self.S):
                # for sto in range(0,self.S, 2):
                    # f_plus = self.V[time-1, sto]
                    # f_minus = self.V[time-1, sto+1]
                    print('sto flag')
                    self.m.Equation(f_plus[0, sto] <= V[0, sto] - self.sto['Vmin'].values[sto])
                    # self.m.Equation(f_minus[0, sto] <= self.sto['Vmax'].values[sto] - V[0, sto])
                    # self.m.Equation(V[time, sto] == (V[time-1, sto] + f_minus[time-1, sto] - f_plus[time-1, sto])) 
            elif self.T >1:
                for time in range(1, self.T):
                    for sto in range(self.S):
                    # for sto in range(0,self.S, 2):
                        # f_plus = self.V[time-1, sto]
                        # f_minus = self.V[time-1, sto+1]
                        self.m.Equation(f_plus[time-1, sto] <= V[time-1, sto] - self.sto['Vmin'].values[sto])
                        self.m.Equation(f_minus[time-1, sto] <= self.sto['Vmax'].values[sto] - V[time-1, sto])
                        self.m.Equation(V[time, sto] == (V[time-1, sto] + f_minus[time-1, sto] - f_plus[time-1, sto])) 
                for sto in range(self.S):
                    self.m.Equation(f_plus[time, sto] <= V[time, sto] - self.sto['Vmin'].values[sto])
                    self.m.Equation(f_minus[time, sto] <= self.sto['Vmax'].values[sto] - V[time, sto])

    def solve_network(self):
        self.m.options.OTOL = 1e-7
        self.m.options.MAX_ITER = 1e9
        self.m.solver_options = ['mu_strategy adaptive',
                                 'constr_viol_tol 1e-7',
                                 'acceptable_tol 1e-7',
                                 'nlp_scaling_method none']
        #  'bound_push 1e-10',\
        #  'bound_frac 1e-10']
        # self.m.Minimize(self.J)
        self.m.open_folder()
        try:
          self.m.solve(disp=self.disp)
        except:
            print('Not successful')
            from gekko.apm import get_file
            print(self.m._server)
            print(self.m._model_name)
            f = get_file(self.m._server,self.m._model_name,'infeasibilities.txt')
            f = f.decode().replace('\r','')
            with open('infeasibilities.txt', 'w') as fl:
                fl.write(str(f))

    def get_values(self, X):
        result = []
        for i in X.flatten():
            try:
                result.append(i.value[0])
            except:
                result.append(i)
        return np.array(result).reshape(X.shape)

    def costs(self, ):
        obj_gas_cost = self.get_values(self.X_gas.flatten()) @ self.cost_gas
        obj_power_cost = self.get_values(self.X_power.flatten()) @ self.cost_power
        return obj_gas_cost[0], obj_power_cost[0], (obj_gas_cost + obj_power_cost)[0]  


#%%
# 118IEEE - 48 node gas 
path_power = '/home/usuario/Work/NGP/db/ps_case9.xlsx'
path_gas = '/home/usuario/Work/NGP/db/ng_case8.xlsx'
name = '9bus_8node'

df_power = pd.ExcelFile(path_power)
df_gas = pd.ExcelFile(path_gas)


savepath = '/home/usuario/Work/NGP/montecarlo_results/'

# types = ['MPCC', 'SOC', 'Taylor']
types = ['Taylor']
samples = 100
np.random.seed(42)
# seeds = np.random.randint(1, 3000,size=(samples))
seeds = [1022, 2614, 1130, 1501, 703, 2450, 2778, 1580, 162, 202, 1982, 
996, 2318, 816, 456, 1276, 1017, 2344, 2768, 338, 879, 1077, 792, 2265, 764,
 2236, 380, 493, 1181, 2063, 65, 2569, 1368, 1153]
T = 1

for typ in types:
        objects = []
        savepath_folder = '/home/usuario/Work/NGP/montecarlo_results/' + typ
        for i in range(samples):
            p = InterconnectedEnergySystemOptimizer(data_gas=df_gas,
                                                    data_power=df_power,
                                                    T=T,
                                                    weymouth_type=typ,
                                                    seed=seeds[i])
            objects.append(p)
            # costs.append(p.costs())
            file_ = open(savepath_folder+'/'+name+typ+'32.pkl', 'wb')
            dill.dump(objects, file_)
            file_.close()
