% Example usage of ryb_colormap function
% This script demonstrates how to use the custom colormap with different step sizes

%% Example 1: Basic usage with default step size
figure(1);
subplot(2,2,1);
data1 = peaks(50);
imagesc(data1);
colormap(ryb_colormap(64));  % Default 64 colors
colorbar;
title('Default (64 colors)');

%% Example 2: High resolution colormap
subplot(2,2,2);
imagesc(data1);
colormap(ryb_colormap(256));  % High resolution
colorbar;
title('High resolution (256 colors)');

%% Example 3: Low resolution colormap
subplot(2,2,3);
imagesc(data1);
colormap(ryb_colormap(16));  % Low resolution
colorbar;
title('Low resolution (16 colors)');

%% Example 4: Medium resolution
subplot(2,2,4);
imagesc(data1);
colormap(customcolormap_preset('red-yellow-blue'));  % Medium-high resolution
colorbar;
title('Medium resolution (128 colors)');

%% Example 5: Visualize the colormap itself
figure(2);
for n = [16, 32, 64, 128]
    subplot(2,2,find([16, 32, 64, 128] == n));
    cmap = ryb_colormap(n);
    image(permute(cmap, [1, 3, 2]));
    axis off;
    title(sprintf('Colormap with %d colors', n));
end

