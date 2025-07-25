import torch.nn.functional as F
from timeit import default_timer
from utilities3 import *
import matplotlib.pyplot as plt
torch.manual_seed(0)
np.random.seed(123)
################################################################
# fourier layer
################################################################

prefix = "Data/"
class SpectralConv2d(nn.Module):
    def __init__(self, in_channels, out_channels, modes1, modes2):
        super(SpectralConv2d, self).__init__()

        """
        2D Fourier layer. It does FFT, linear transform, and Inverse FFT.
        """

        self.in_channels = in_channels
        self.out_channels = out_channels
        self.modes1 = modes1 #Number of Fourier modes to multiply, at most floor(N/2) + 1
        self.modes2 = modes2

        self.scale = (1 / (in_channels * out_channels))
        self.weights1 = nn.Parameter(self.scale * torch.rand(in_channels, out_channels, self.modes1, self.modes2, dtype=torch.cfloat))
        self.weights2 = nn.Parameter(self.scale * torch.rand(in_channels, out_channels, self.modes1, self.modes2, dtype=torch.cfloat))

    # Complex multiplication
    def compl_mul2d(self, input, weights):
        # (batch, in_channel, x,y ), (in_channel, out_channel, x,y) -> (batch, out_channel, x,y)
        return torch.einsum("bixy,ioxy->boxy", input, weights)

    def forward(self, x):
        batchsize = x.shape[0]
        #Compute Fourier coeffcients up to factor of e^(- something constant)
        x_ft = torch.fft.rfft2(x)

        # Multiply relevant Fourier modes
        out_ft = torch.zeros(batchsize, self.out_channels,  x.size(-2), x.size(-1)//2 + 1, dtype=torch.cfloat, device=x.device)
        out_ft[:, :, :self.modes1, :self.modes2] = \
            self.compl_mul2d(x_ft[:, :, :self.modes1, :self.modes2], self.weights1)
        out_ft[:, :, -self.modes1:, :self.modes2] = \
            self.compl_mul2d(x_ft[:, :, -self.modes1:, :self.modes2], self.weights2)

        #Return to physical space
        x = torch.fft.irfft2(out_ft, s=(x.size(-2), x.size(-1)))
        return x

class MLP(nn.Module):
    def __init__(self, in_channels, out_channels, mid_channels):
        super(MLP, self).__init__()
        self.mlp1 = nn.Conv2d(in_channels, mid_channels, 1)
        self.mlp2 = nn.Conv2d(mid_channels, out_channels, 1)

    def forward(self, x):
        x = self.mlp1(x)
        x = F.gelu(x)
        x = self.mlp2(x)
        return x

class FNO2d(nn.Module):
    def __init__(self, modes1, modes2,  width):
        super(FNO2d, self).__init__()

        """
        The overall network. It contains 4 layers of the Fourier layer.
        1. Lift the input to the desire channel dimension by self.fc0 .
        2. 4 layers of the integral operators u' = (W + K)(u).
            W defined by self.w; K defined by self.conv .
        3. Project from the channel space to the output space by self.fc1 and self.fc2 .

        input: the solution of the coefficient function and locations (a(x, y), x, y)
        input shape: (batchsize, x=s, y=s, c=3)
        output: the solution
        output shape: (batchsize, x=s, y=s, c=1)
        """

        self.modes1 = modes1
        self.modes2 = modes2
        self.width = width
        self.padding = 9 # pad the domain if input is non-periodic

        self.p = nn.Linear(3, self.width) # input channel is 3: (a(x, y), x, y)
        self.conv0 = SpectralConv2d(self.width, self.width, self.modes1, self.modes2)
        self.conv1 = SpectralConv2d(self.width, self.width, self.modes1, self.modes2)
        self.conv2 = SpectralConv2d(self.width, self.width, self.modes1, self.modes2)
        self.conv3 = SpectralConv2d(self.width, self.width, self.modes1, self.modes2)
        self.mlp0 = MLP(self.width, self.width, self.width)
        self.mlp1 = MLP(self.width, self.width, self.width)
        self.mlp2 = MLP(self.width, self.width, self.width)
        self.mlp3 = MLP(self.width, self.width, self.width)
        self.w0 = nn.Conv2d(self.width, self.width, 1)
        self.w1 = nn.Conv2d(self.width, self.width, 1)
        self.w2 = nn.Conv2d(self.width, self.width, 1)
        self.w3 = nn.Conv2d(self.width, self.width, 1)
        self.q = MLP(self.width, 1, self.width * 4) # output channel is 1: u(x, y)

    def forward(self, x):
        grid = self.get_grid(x.shape, x.device)
        x = torch.cat((x, grid), dim=-1)
        x = self.p(x)
        x = x.permute(0, 3, 1, 2)
        x = F.pad(x, [0,self.padding, 0,self.padding])

        x1 = self.conv0(x)
        x1 = self.mlp0(x1)
        x2 = self.w0(x)
        x = x1 + x2
        x = F.gelu(x)

        x1 = self.conv1(x)
        x1 = self.mlp1(x1)
        x2 = self.w1(x)
        x = x1 + x2
        x = F.gelu(x)

        x1 = self.conv2(x)
        x1 = self.mlp2(x1)
        x2 = self.w2(x)
        x = x1 + x2
        x = F.gelu(x)

        x1 = self.conv3(x)
        x1 = self.mlp3(x1)
        x2 = self.w3(x)
        x = x1 + x2

        x = x[..., :-self.padding, :-self.padding]
        x = self.q(x)
        x = x.permute(0, 2, 3, 1)
        return x

    def get_grid(self, shape, device):
        batchsize, size_x, size_y = shape[0], shape[1], shape[2]
        gridx = torch.tensor(np.linspace(0, 1, size_x), dtype=torch.float)
        gridx = gridx.reshape(1, size_x, 1, 1).repeat([batchsize, 1, size_y, 1])
        gridy = torch.tensor(np.linspace(0, 1, size_y), dtype=torch.float)
        gridy = gridy.reshape(1, 1, size_y, 1).repeat([batchsize, size_x, 1, 1])
        return torch.cat((gridx, gridy), dim=-1).to(device)

ntest, batch_size, lr, epochs = 8, 1, 1e-3, 100
modes = width = 16
s1, s2 = 140, 100

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

a0 = np.loadtxt(prefix + 'vpd_matrix.txt').reshape(s1, s2, 18).transpose(2, 0, 1)
aT = a0.copy()

def run_once(ntrain: int) -> float:
    torch.manual_seed(0)
    np.random.seed(123)

    x_tr = torch.from_numpy(a0[:ntrain].reshape(-1).astype(np.float32)).to(device)
    y_tr = torch.from_numpy(aT[1:ntrain+1].reshape(-1).astype(np.float32)).to(device)
    x_te = torch.from_numpy(a0[ntrain:ntrain+ntest].reshape(-1).astype(np.float32)).to(device)
    y_te = torch.from_numpy(aT[ntrain+1:ntrain+ntest+1].reshape(-1).astype(np.float32)).to(device)

    x_tr = x_tr.view(ntrain, s1, s2, 1)
    y_tr = y_tr.view(ntrain, s1, s2, 1)
    x_te = x_te.view(ntest , s1, s2, 1)
    y_te = y_te.view(ntest , s1, s2, 1)

    train_loader = torch.utils.data.DataLoader(
        torch.utils.data.TensorDataset(x_tr, y_tr), batch_size, shuffle=True)
    test_loader = torch.utils.data.DataLoader(
        torch.utils.data.TensorDataset(x_te, y_te), batch_size, shuffle=False)

    model = FNO2d(modes, modes, width).to(device)
    opt = torch.optim.Adam(model.parameters(), lr, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(
        opt, T_max=epochs * (ntrain // batch_size))
    loss_fn = LpLoss(size_average=False)

    for _ in range(epochs):
        model.train()
        for x, y in train_loader:
            opt.zero_grad()
            out = model(x).reshape(batch_size, s1, s2)
            loss = loss_fn(out.view(batch_size, -1), y.view(batch_size, -1))
            loss.backward()
            opt.step()
            sched.step()

        model.eval()
        with torch.no_grad():
            for x, y in test_loader:
                _ = model(x)

    test_l2 = 0.0
    model.eval()
    with torch.no_grad():
        for x, y in test_loader:
            out = model(x).reshape(batch_size, s1, s2)
            test_l2 += ((out.view(1, -1) - y.view(1, -1)) ** 2).sum().item()

    return (test_l2 / (s1 * s2) / ntest) ** 0.5

results = []
for ntrain in range(4, 10):
    rmse = run_once(ntrain)
    results.append(rmse)
    print(f'ntrain={ntrain}, RMSE={rmse:.6f}')

print('All RMSE:', results)

np.savetxt(prefix+"rolling_window_FNO.txt", results, fmt="%.6f")