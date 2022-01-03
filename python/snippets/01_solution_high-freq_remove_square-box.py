# 0 entries: points not sampled, 1 entries: points to be sampled
# Initialize k-space with "1" everywhere and then place the "0" appropriately
kspace_maskbox = np.ones((img_size,img_size), dtype="float64")
# use fancy indexing along rows
kspace_maskbox[idx_vec_, ] = 0

list_img_size = np.arange(0., img_size).tolist()
filtered_center = [x for x in list_img_size if x not in idx_vec_]
array_idx_center = np.array(filtered_center)
array_idx_center_ = array_idx_center.astype("int")
# use fancy indexing along cols
kspace_maskbox[:, array_idx_center_] = 1

# Note that the combination of these two fancy indexing replaces the poor code below
#for i in idx_vec_:
#    for j in idx_vec_:
#        kspace_mask[i,j]=0.

# Generate the kspace data: first Fourier transform the image
kspace_data = np.fft.fftshift(fft(mri_img))
#add Gaussian complex-valued random noise
signoise = 10
kspace_data += np.random.randn(*mri_img.shape) * signoise * (1+1j)
# Mask data to perform subsampling
kspace_data *= kspace_maskbox 

# Zero order solution
image_rec0 = ifft(np.fft.ifftshift(kspace_data))

fig, axs = plt.subplots(2, 2, figsize=(10, 10) )
axs[0,0].imshow(mri_img, cmap='Greys_r')
axs[0,0].set_title("True image")
axs[0,1].imshow(kspace_maskbox, cmap='Greys_r')
axs[0,1].set_title("High frequency box sampling mask")
axs[1,0].imshow(np.abs(kspace_data),  cmap='gray', vmax=1*np.abs(kspace_data).max())
#axs[1].imshow(np.abs(np.fft.ifftshift(kspace_data)), cmap='Greys_r')
axs[1,0].set_title("k-space noisy data")
axs[1,1].imshow(np.abs(image_rec0), cmap='Greys_r')
axs[1,1].set_title("Zero-order recon")
plt.show()

