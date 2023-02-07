#code cloned from https://github.com/eriklindernoren/Keras-GAN/blob/master/cyclegan

#modified for thermal-to-visual application

import scipy
from glob import glob
import numpy as np

aSet="thermal_bmp"
#bSet="classic"
bSet="infrared"
class DataLoader():
    def __init__(self, dataset_name, img_res=(128, 128)):
        self.dataset_name = dataset_name
        self.img_res = img_res

    def load_data(self, domain, dataName, batch_size=1, is_testing=False):
        if is_testing:
            if domain=="A":
                path = glob('carl_database/%s/*/session_4/%s/*' % (dataName,aSet))
            else:
                path = glob('carl_database/%s/*/session_4/%s/*' % (dataName,bSet))
        else:
            if domain=="A":
                path = glob('carl_database/%s/*/session_[1-4]/%s/*' % (self.dataset_name,aSet))
            else:
                path = glob('carl_database/%s/*/session_[1-4]/%s/*' % (self.dataset_name,bSet))

        if batch_size==0:
            batch_images = path
        else:
            if is_testing==True:
                batch_images = path[0::int(len(path)/batch_size)]
                batch_images = batch_images[0:batch_size]
            else:
                batch_images = np.random.choice(path, size=batch_size)

        imgs = []
        for img_path in batch_images:
            if(img_path[-2:]=='db'):
                continue
            img = self.imread(img_path)

            if not is_testing:
                img = scipy.misc.imresize(img, self.img_res)
            else:
                img = scipy.misc.imresize(img, self.img_res)

            imgs.append(img)

        #normalize pixel values to be between -1 and 1
        imgs = np.array(imgs)/127.5 - 1.

        return imgs

    def load_batch(self, batch_size=1, is_testing=False):
        #load from training or testing set
        if is_testing:
            path_A = glob('carl_database/%s/*/session_4/%s/*' % (self.dataset_name,aSet))
            path_B = glob('carl_database/%s/*/session_4/%s/*' % (self.dataset_name,bSet))
        else:
            path_A = glob('carl_database/%s/*/session_[1-4]/%s/*' % (self.dataset_name,aSet))
            path_B = glob('carl_database/%s/*/session_[1-4]/%s/*' % (self.dataset_name,bSet))


        self.n_batches = int(min(len(path_A), len(path_B)) / batch_size)
        total_samples = self.n_batches * batch_size

        path_A = np.random.choice(path_A, total_samples, replace=False)
        path_B = np.random.choice(path_B, total_samples, replace=False)

        for i in range(self.n_batches-1):
            batch_A = path_A[i*batch_size:(i+1)*batch_size]
            batch_B = path_B[i*batch_size:(i+1)*batch_size]
            imgs_A, imgs_B = [], []
            for img_A, img_B in zip(batch_A, batch_B):
                if(img_A[-2:]=='db' or img_B[-2:]=='db'):
                    continue
                img_A = self.imread(img_A)
                img_B = self.imread(img_B)

                img_A = scipy.misc.imresize(img_A, self.img_res)
                img_B = scipy.misc.imresize(img_B, self.img_res)

                imgs_A.append(img_A)
                imgs_B.append(img_B)

            imgs_A = np.array(imgs_A)/127.5 - 1.
            imgs_B = np.array(imgs_B)/127.5 - 1.

            yield imgs_A, imgs_B


    def imread(self, path):

        import cv2
        im = cv2.imread(path).astype(np.float)
        return im