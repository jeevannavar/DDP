import numpy as np
import torch
import torch.nn.functional as F

cuda = True if torch.cuda.is_available() else False


def cal_sample_weight(labels, num_class, use_sample_weight=True):
    """ calculate sample weights based on training label distribution
    """
    if not use_sample_weight:
        return np.ones(len(labels)) / len(labels)
    count = np.zeros(num_class)
    for i in range(num_class):
        count[i] = np.sum(labels==i)
    sample_weight = np.zeros(labels.shape)
    for i in range(num_class):
        sample_weight[np.where(labels==i)[0]] = count[i]/np.sum(count)
    
    return sample_weight


def one_hot_tensor(y, num_dim):
    """convet y to one hot coding"""
    y_onehot = torch.zeros(y.shape[0], num_dim)
    y_onehot.scatter_(1, y.view(-1,1), 1)
    
    return y_onehot


def cosine_distance_torch(x1, x2=None, eps=1e-8):
    """cosine distance with torch"""
    x2 = x1 if x2 is None else x2
    w1 = x1.norm(p=2, dim=1, keepdim=True)
    w2 = w1 if x2 is x1 else x2.norm(p=2, dim=1, keepdim=True)
    return 1 - torch.mm(x1, x2.t()) / (w1 * w2.t()).clamp(min=eps)


def to_sparse(x):
    """ converts dense tensor x to sparse format """
    x_typename = torch.typename(x).split('.')[-1]
    sparse_tensortype = getattr(torch.sparse, x_typename)

    indices = torch.nonzero(x)
    if len(indices.shape) == 0:  # if all elements are zeros
        return sparse_tensortype(*x.shape)
    indices = indices.t()
    values = x[tuple(indices[i] for i in range(indices.shape[0]))]
    return sparse_tensortype(indices, values, x.size())


def cal_adj_mat_parameter(edge_per_node, data, metric="cosine"):
    """ generating parameter for the weighted adjacency matrix adaptively
        according to average edge per node (does not include self loop)
        for training data
    """
    assert metric == "cosine", "Only cosine distance implemented"
    dist = cosine_distance_torch(data, data)
    parameter = torch.sort(dist.reshape(-1,)).values[edge_per_node*data.shape[0]]
    return np.asscalar(parameter.data.cpu().numpy())


def graph_from_dist_tensor(dist, parameter, self_dist=True):
    """
    generate graph edge mask from distance matrix
    parameter: radius threshold
    self_dist: is dist self-pairwise distance
    Return: g
    """
    if self_dist:
        assert dist.shape[0]==dist.shape[1], "Input is not pairwise dist matrix"
    # retrain edge with distance smaller than parameter
    g = (dist <= parameter).float()
    if self_dist: # square adj mat
        diag_idx = np.diag_indices(g.shape[0])
        g[diag_idx[0], diag_idx[1]] = 0
        
    return g


def gen_adj_mat_tensor(data, parameter, metric="cosine"):
    """
    generate weighted graph adjacency matrix for GCN training
    parameter: radius threshold
    Return: adj (Tensor)
    """
    assert metric == "cosine", "Only cosine distance implemented"
    # calculate pairwise distance
    dist = cosine_distance_torch(data, data)
    # construct binary graph
    g = graph_from_dist_tensor(dist, parameter, self_dist=True)

    # convert to affinity matrix
    if metric == "cosine":
        adj = 1-dist
    else:
        raise NotImplementedError
    adj = adj*g # retain selected edges

    adj_T = adj.transpose(0,1)
    I = torch.eye(adj.shape[0])
    if cuda:
        I = I.cuda()
    # build symmetric adjacency matrix
    adj = adj + adj_T*(adj_T > adj).float() - adj*(adj_T > adj).float()
    # add self-cycle, row normalization
    adj = F.normalize(adj + I, p=1)
    # to sparse COOrdinate format
    adj = to_sparse(adj)
    
    return adj


def gen_test_adj_mat_tensor(data, trte_idx, parameter, metric="cosine"):
    """ generate weighted graph adjacency matrix for testing
    """
    assert metric == "cosine", "Only cosine distance implemented"
    adj = torch.zeros((data.shape[0], data.shape[0]))
    if cuda:
        adj = adj.cuda()
    num_tr = len(trte_idx["tr"])
    
    # Training data to testing data
    dist_tr2te = cosine_distance_torch(data[trte_idx["tr"]], data[trte_idx["te"]])
    g_tr2te = graph_from_dist_tensor(dist_tr2te, parameter, self_dist=False)
    if metric == "cosine":
        adj[:num_tr,num_tr:] = 1-dist_tr2te
    else:
        raise NotImplementedError
    adj[:num_tr,num_tr:] = adj[:num_tr,num_tr:]*g_tr2te
    
    # Testing data to training data
    dist_te2tr = cosine_distance_torch(data[trte_idx["te"]], data[trte_idx["tr"]])
    g_te2tr = graph_from_dist_tensor(dist_te2tr, parameter, self_dist=False)
    if metric == "cosine":
        adj[num_tr:,:num_tr] = 1-dist_te2tr
    else:
        raise NotImplementedError
    adj[num_tr:,:num_tr] = adj[num_tr:,:num_tr]*g_te2tr # retain selected edges
    
    adj_T = adj.transpose(0,1)
    I = torch.eye(adj.shape[0])
    if cuda:
        I = I.cuda()
    # build symmetric adjacency matrix
    adj = adj + adj_T*(adj_T > adj).float() - adj*(adj_T > adj).float()
    # add self-cycle, row normalization
    adj = F.normalize(adj + I, p=1)
    # to sparse COOrdinate format
    adj = to_sparse(adj)
    
    return adj


def print_dict(dictionary, print_format=".3f"):
    """ generate string based on dictionary content key:value"""
    output_list = []
    for key in dictionary:
        output_list.append("{:}:{:{:}}".format(key, dictionary[key], print_format))
    return ", ".join(output_list)