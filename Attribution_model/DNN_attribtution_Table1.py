import os
import torch.nn.functional as F
import torch
from timeit import default_timer
import matplotlib.pyplot as plt
torch.manual_seed(0)

import numpy as np
import torch.nn as nn
import torch.optim as optim

np.random.seed(123)

class RNNModel(nn.Module):
    def __init__(self):
        super(RNNModel, self).__init__()
        self.rnn = nn.RNN(input_size=2, hidden_size=64, batch_first=True)
        self.fc1 = nn.Linear(64, 64)  
        self.fc2 = nn.Linear(64, 32)
        self.fc3 = nn.Linear(32, 1)   
    def forward(self, x):
        x, _ = self.rnn(x)  
        x = torch.relu(self.fc1(x)) 
        x = torch.relu(self.fc2(x))
        x = self.fc3(x)  
        return x



class DNNModel(nn.Module):
    def __init__(self):
        super(DNNModel, self).__init__()
        # Assuming the input size to the network is 2
        self.fc1 = nn.Linear(2, 64)
        self.fc2 = nn.Linear(64, 64)
        self.fc3 = nn.Linear(64, 1)

    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        x = self.fc3(x)
        return x




prefix = os.path.join("Data")


s=14000


# prepare the data

vpd = np.loadtxt(os.path.join(prefix, 'vpd_matrix.txt'))
precip = np.loadtxt(os.path.join(prefix, 'precip_matrix.txt'))
ndvi = np.loadtxt(os.path.join(prefix, 'ndvi_matrix.txt'))


vpd=vpd.reshape(s, 18)
precip=precip.reshape(s, 18)
ndvi=ndvi.reshape(s, 18)




vpd = torch.tensor(vpd, dtype=torch.float32).unsqueeze(2) 
precip = torch.tensor(precip, dtype=torch.float32).unsqueeze(2) 
ndvi = torch.tensor(ndvi, dtype=torch.float32).unsqueeze(2)  

features = torch.cat((vpd, precip), dim=2)  

device = 'cuda' if torch.cuda.is_available() else 'cpu'
all_predictions= torch.zeros(s, 8, dtype=torch.float32).to(device)


l2_losses = []

for i in range(14000):
    model = DNNModel().to(device)
    optimizer = optim.Adam(model.parameters(), lr=0.001)
    criterion = nn.MSELoss()

    grid_features = features[i].view(1, -1, 2).to(device) 
    grid_ndvi = ndvi[i].view(1, -1, 1).to(device)  
    
    train_data = grid_features[:, :10, :]
    test_data = grid_features[:, 10:, :]
    train_targets = grid_ndvi[:, :10, :]
    test_targets = grid_ndvi[:, 10:, :]


    for epoch in range(50): 
        model.train()
        optimizer.zero_grad()
        outputs = model(train_data)
        loss = criterion(outputs, train_targets)
        loss.backward()
        optimizer.step()

    model.eval()
    with torch.no_grad():
        test_pred = model(test_data)
        all_predictions[i] = test_pred.view(-1) 
        mse_criterion = torch.nn.MSELoss()
        l2_loss = mse_criterion(test_pred, test_targets)
        l2_losses.append(l2_loss.item()) 
        print(f'Grid {i} Test Loss (L2/MSE): {l2_loss.item()}')



print((sum(l2_losses)/14000)**0.5)
