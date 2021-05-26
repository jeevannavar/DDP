""" Componets of the model
"""
import torch.nn as nn
import torch
import torch.nn.functional as F


def xavier_init(m):
    """
    Xavier weight initialization for fully-connected layer
    """
    if type(m) == nn.Linear:
        nn.init.xavier_normal_(m.weight)
        if m.bias is not None:
           m.bias.data.fill_(0.0)
           

class GraphConvolution(nn.Module):
    """#===========================================
    Simple GCN layer
    https://github.com/meliketoy/graph-cnn.pytorch/blob/master/layers.py
    """#===========================================

    def __init__(self, in_features, out_features, bias=True):
        super().__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.weight = nn.Parameter(torch.FloatTensor(in_features, out_features))
        if bias:
            self.bias = nn.Parameter(torch.FloatTensor(out_features))
        # initialization
        nn.init.xavier_normal_(self.weight.data)
        if self.bias is not None:
            self.bias.data.fill_(0.0)
    
    def forward(self, x, adj):
        support = torch.mm(x, self.weight)
        output = torch.sparse.mm(adj, support)
        if self.bias is not None:
            return output + self.bias
        else:
            return output
    

class GCN_E_2(nn.Module):
    """#===========================================
        GCN encoder with two layers
    """#===========================================
    
    def __init__(self, in_dim, hgcn_dim, out_dim, dropout):
        super().__init__()
        self.gc1 = GraphConvolution(in_dim, hgcn_dim)
        self.gc2 = GraphConvolution(hgcn_dim, out_dim)
        self.dropout = dropout

    def forward(self, x, adj):
        x = self.gc1(x, adj)
        x = F.leaky_relu(x, 0.25)
        x = F.dropout(x, self.dropout, training=self.training)
        x = self.gc2(x, adj)
        x = F.leaky_relu(x, 0.25)
        x = F.dropout(x, self.dropout, training=self.training)

        return x

class GCN_E_N(nn.Module):
    """#===========================================
        GCN encoder with two layers
    """#===========================================
    
    def __init__(self, in_dim, hgcn_dim, out_dim, dropout, num_layers):
        super().__init__()
        self.gc = [GraphConvolution(in_dim, hgcn_dim)]
        for _ in range(num_layers-2):
            self.gc.append(GraphConvolution(hgcn_dim, hgcn_dim))
        self.gc.append(GraphConvolution(hgcn_dim, out_dim))
        self.dropout = dropout

    def forward(self, x, adj):
        for i in range(len(self.gc)):
            x = self.gc[i](x, adj)
            x = F.leaky_relu(x, 0.25)
            x = F.dropout(x, self.dropout, training=self.training)

        return x

class Classifier_1(nn.Module):
    """#===========================================
        One layer fully connected classifier
    """#===========================================

    def __init__(self, in_dim, out_dim):
        super().__init__()
        self.clf = nn.Sequential(nn.Linear(in_dim, out_dim))
        self.clf.apply(xavier_init)

    def forward(self, x):
        x = self.clf(x)
        return x


class VCDN(nn.Module):
    """#===========================================
        View-correlation discovery network from "Generative Multi-View Human Action Recognition"
        Input: predicted logits from GCN classifiers
        Output: final predicted logits
    """#===========================================

    def __init__(self, num_cls, hvcdn_dim, num_views):
        """
        Initialization:
        Inputs:
            num_cls -> number of classes/logits
            hvcdn_dim -> input dimension for the vcdn, which is num_cls**num_views
            num_views -> number of views/gcns
        """
        super().__init__()
        self.num_cls = num_cls
        self.model = nn.Sequential(
            nn.Linear(num_cls**num_views, hvcdn_dim),
            nn.LeakyReLU(0.25),
            nn.Linear(hvcdn_dim, num_cls)
        )
        self.model.apply(xavier_init)
        
    def forward(self, views):
        """
        Inputs:
            views -> a list of gcns whose output is being supplied to the vcdn
        Outputs:
            output -> labels
        """
        for i in range(len(views)):
            views[i] = torch.sigmoid(views[i])
        if len(views) == 1:
            output = self.model(views[0])
        else:
            x1 = views[0].unsqueeze(-1)
            x2 = views[1].unsqueeze(1)
            vcdn_mat = torch.reshape(torch.matmul(x1, x2),(-1,self.num_cls*self.num_cls,1))

            for i in range(2, len(views)):
                x_n = views[i].unsqueeze(1)
                vcdn_mat = torch.reshape(torch.matmul(vcdn_mat, x_n),(-1,self.num_cls**(i+1),1))
            vcdn_feat = torch.reshape(vcdn_mat, (-1,self.num_cls**len(views)))
            output = self.model(vcdn_feat)

        return output

class DenseCombine(nn.Module):
    """#===========================================
        A pair of fully-connected layers to combine the outputs from individual GCNs
        Input: predicted logits from GCN classifiers
        Output: final predicted logits
    """#===========================================

    def __init__(self, num_cls, num_views):
        """
        Initialization:
        Inputs:
            num_cls -> number of classes/logits
            hvcdn_dim -> input dimension for the vcdn, which is num_cls**num_views
            num_views -> number of views/gcns
        """
        super().__init__()
        self.num_cls = num_cls
        self.model = nn.Sequential(
            nn.Linear(num_cls*num_views, num_cls*num_views),
            nn.LeakyReLU(0.25),
            nn.Linear(num_cls*num_views, num_cls)
        )
        self.model.apply(xavier_init)
        
    def forward(self, views):
        """
        Inputs:
            views -> a list of individual gcn outputs
        Outputs:
            output -> labels
        """

        for i in range(len(views)):
            views[i] = torch.sigmoid(views[i])

        else:
            concat = torch.cat(views, dim=1)
            output = self.model(concat)

        return output
    
def init_model_dict(num_class, dim_list, dim_he_list, dim_hc, gcn_names, combiner, gcn_dropout=0.5, num_gcn=2):
    """
        Initialize the models - GCNs and VCDN
        Inputs:
            num_class   -> number of classes/logits
            dim_list    -> list of input dimensions for each GCN
            dim_he_list -> list of hidden layer dimensions for each GCN
            dim_hc      -> input dimension for the VCDN
            gcn_names   -> list of namesfor the gcns and the combiner, if it is being used
            combiner    -> name of combiner, false if no combiner is used
            gcn_dropout -> dropout rate for the GCN layers
            num_gcn     -> number of gcn layers to train on

        Outputs:
            model_dict  -> dictionary of individual gcn and vcdn layers/models
                            It contains keys in the format:
                                "Ei" for GCN encoder layer i
                                gcn_name for GCN classifier layer
                                combiner for combiner
    """
    model_dict = {}
    for i in range(len(dim_list)):
        encoder = f"E{i+1}"
        model_dict[encoder] = GCN_E_N(dim_list[i], dim_he_list[i], dim_he_list[i], gcn_dropout, num_gcn)
        classifier = gcn_names[i]
        model_dict[classifier] = Classifier_1(dim_he_list[i], num_class)
    
    num_views = len(dim_list)

    if combiner:
        if combiner == "VCDN":
            model_dict[combiner] = VCDN(num_class, dim_hc, num_views)
        if combiner == "FullyConnected":
            model_dict[combiner] = DenseCombine(num_class, num_views)
    
    return model_dict


def init_optim(model_dict, dim_list, gcn_names, combiner, lr=1e-3, weight_decay=1e-4):
    """
        Initialize the optimization dictionary
        Inputs:
            model_dict -> dictionary containing all the individual GCNs and the VCDN
            lr         -> learning rate
            dim_list   -> list of input dimensions for the GCNs
            gcn_names  -> list of names for the gcns and the combiner, if it is being used
            combiner   -> name of combiner, false if no combiner is used
        Outputs:
            optim_dict -> Adam optimiser for each of the GCN and VCDN layers
    """
    optim_dict = {}
    for i in range(len(dim_list)):
        classifier = gcn_names[i]
        encoder = f"E{i+1}"
        optim_dict[classifier] = torch.optim.Adam(list(model_dict[encoder].parameters())+list(model_dict[classifier].parameters()), lr=lr, weight_decay=weight_decay)

    if combiner:
        optim_dict[combiner] = torch.optim.Adam(model_dict[combiner].parameters(), lr=lr, weight_decay=weight_decay)
    
    return optim_dict