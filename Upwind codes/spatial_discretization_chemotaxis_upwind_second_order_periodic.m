function [vu_x, vu_y] = spatial_discretization_chemotaxis_upwind_second_order_periodic(u, v, hx, hy, nx, ny)
% Computes the spatial discretization for the chemotaxis term with second-order
% upwind schemes and periodic boundary conditions.
%
% Inputs:
% u - cell density matrix (nx x ny)
% v - chemical concentration matrix (nx x ny)
% hx, hy - grid spacing in the x and y directions
% nx, ny - number of grid points in the x and y directions
%
% Outputs:
% vu_x, vu_y - discretized chemotaxis terms in the x and y directions

% Initialize chemotaxis terms
vu_x = zeros(nx, ny);
vu_y = zeros(nx, ny);
v_x = zeros(nx, ny);
v_y = zeros(nx, ny);

% Compute central differences for v gradients (v_x and v_y)

for i = 1:nx
    for j = 1:ny
        % Calculate v_x with periodic boundary conditions in the x direction
        if i == 1
            v_x(1, j) = (v(2, j) - v(nx, j)) / (2 * hx);
        elseif i == nx
            v_x(nx, j) = (v(1, j) - v(nx-1, j)) / (2 * hx);
        else
            v_x(i, j) = (v(i+1, j) - v(i-1, j)) / (2 * hx);
        end

        % Calculate v_y with periodic boundary conditions in the y direction
        if j == 1
            v_y(i, 1) = (v(i, 2) - v(i, ny)) / (2 * hy);
        elseif j == ny
            v_y(i, ny) = (v(i, 1) - v(i, ny-1)) / (2 * hy);
        else
            v_y(i, j) = (v(i, j+1) - v(i, j-1)) / (2 * hy);
        end
    end
end

% Compute chemotaxis terms with second-order upwind schemes
for i = 1:nx
    for j = 1:ny
        % Chemotaxis term in the x-direction
        if v_x(i, j) >= 0
            % v_x >= 0: Second-order upwind scheme in x-direction
            if i > 2
                vu_x(i, j) = (1/hx) * ((3/2) * v_x(i, j) * u(i, j) - 2 * v_x(i-1, j) * u(i-1, j) + (1/2) * v_x(i-2, j) * u(i-2, j));
            else
                % Periodic boundary adjustment for i = 1, 2
                if i == 2
                    vu_x(i, j) = (1/hx) * ((3/2) * v_x(i, j) * u(i, j) - 2 * v_x(i-1, j) * u(i-1, j) + (1/2) * v_x(nx, j) * u(nx, j));
                else % i == 1
                    vu_x(i, j) = (1/hx) * (v_x(i, j) * u(i, j) - v_x(nx, j) * u(nx, j));
                end
            end
        else
            % v_x < 0: Second-order upwind scheme in x-direction
            if i < nx-1
                vu_x(i, j) = (1/hx) * (-(1/2) * v_x(i+2, j) * u(i+2, j) + 2 * v_x(i+1, j) * u(i+1, j) - (3/2) * v_x(i, j) * u(i, j));
            else
                % Periodic boundary adjustment for i = nx-1, nx
                if i == nx-1
                    vu_x(i, j) = (1/hx) * (-(1/2) * v_x(1, j) * u(1, j) + 2 * v_x(i+1, j) * u(i+1, j) - (3/2) * v_x(i, j) * u(i, j));
                else % i == nx
                    vu_x(i, j) = (1/hx) * (v_x(1, j) * u(1, j) - v_x(i, j) * u(i, j));
                end
            end
        end

        % Chemotaxis term in the y-direction
        if v_y(i, j) >= 0
            % v_y >= 0: Second-order upwind scheme in y-direction
            if j > 2
                vu_y(i, j) = (1/hy) * ((3/2) * v_y(i, j) * u(i, j) - 2 * v_y(i, j-1) * u(i, j-1) + (1/2) * v_y(i, j-2) * u(i, j-2));
            else
                % Periodic boundary adjustment for j = 1, 2
                if j == 2
                    vu_y(i, j) = (1/hy) * ((3/2) * v_y(i, j) * u(i, j) - 2 * v_y(i, j-1) * u(i, j-1) + (1/2) * v_y(i, ny) * u(i, ny));
                else % j == 1
                    vu_y(i, j) = (1/hy) * (v_y(i, j) * u(i, j) - v_y(i, ny) * u(i, ny));
                end
            end
        else
            % v_y < 0: Second-order upwind scheme in y-direction
            if j < ny-1
                vu_y(i, j) = (1/hy) * (-(1/2) * v_y(i, j+2) * u(i, j+2) + 2 * v_y(i, j+1) * u(i, j+1) - (3/2) * v_y(i, j) * u(i, j));
            else
                % Periodic boundary adjustment for j = ny-1, ny
                if j == ny-1
                    vu_y(i, j) = (1/hy) * (-(1/2) * v_y(i, 1) * u(i, 1) + 2 * v_y(i, j+1) * u(i, j+1) - (3/2) * v_y(i, j) * u(i, j));
                else % j == ny
                    vu_y(i, j) = (1/hy) * (v_y(i, 1) * u(i, 1) - v_y(i, j) * u(i, j));
                end
            end
        end
    end
end

% Flatten output for consistency
vu_x = vu_x(:);
vu_y = vu_y(:);
end
