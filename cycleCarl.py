#code cloned from https://github.com/eriklindernoren/Keras-GAN/blob/master/cyclegan

#modified for thermal-to-visual application

from __future__ import print_function, division
import scipy

from keras.datasets import mnist
from instance_normalization import InstanceNormalization
from keras.layers import Input, Dense, Reshape, Flatten, Dropout, Concatenate
from keras.layers import BatchNormalization, Activation, ZeroPadding2D, InputSpec, Layer
from keras.layers.advanced_activations import LeakyReLU
from keras.layers.convolutional import UpSampling2D, Conv2D
from keras.models import Sequential, Model
from keras.optimizers import Adam
import datetime
import matplotlib.pyplot as plt
import sys
from data_loader2 import DataLoader
import numpy as np
import os
from keras import backend as K
import tensorflow as tf
from keras_contrib.losses import DSSIMObjective

class CycleGAN():
    def __init__(self):
        # Input shape
        self.img_rows = 128
        self.img_cols = 128
        self.channels = 3
        self.img_shape = (self.img_rows, self.img_cols, self.channels)

        self.dataset_name = 'carl_database'
        self.data_loader = DataLoader(dataset_name=self.dataset_name,
                                      img_res=(self.img_rows, self.img_cols))

        self.disc_patch = (8, 8, 1)

        # Number of filters in the first layer of G and D
        self.gf = 32
        self.df = 64

        # Loss weights
        self.lambda_cycle = 10.0                   # Cycle-consistency loss
        self.lambda_id = 0.1 * self.lambda_cycle    # Identity loss

        optimizer = Adam(0.0002, 0.5)

        # Build and compile the discriminators
        self.d_A = self.build_discriminator()
        self.d_B = self.build_discriminator()
        def lse(y_true, y_pred):
            loss = tf.reduce_mean(tf.squared_difference(y_pred, y_true))
            return loss

        def cycle_loss(y_true, y_pred):
            loss = tf.reduce_mean(tf.abs(y_pred - y_true))
            return loss
        loss_weights_D=[0.5]
        self.d_A.compile(loss=['mse'],
            optimizer=optimizer,
            metrics=['accuracy'],loss_weights=loss_weights_D)
        self.d_B.compile(loss=['mse'],
            optimizer=optimizer,
            metrics=['accuracy'],loss_weights=loss_weights_D)

        #-------------------------
        # Construct Computational
        #   Graph of Generators
        #-------------------------

        # Build the generators
        self.g_AB = self.build_generator()
        self.g_BA = self.build_generator()

        # Input images from both domains
        img_A = Input(shape=self.img_shape)
        img_B = Input(shape=self.img_shape)

        # Translate images to the other domain
        fake_B = self.g_AB(img_A)
        fake_A = self.g_BA(img_B)

        # Translate images back to original domain
        reconstr_A = self.g_BA(fake_B)
        reconstr_B = self.g_AB(fake_A)

        # Identity mapping of images
        img_A_id = self.g_BA(img_A)
        img_B_id = self.g_AB(img_B)

        # For the combined model we will only train the generators
        self.d_A.trainable = False
        self.d_B.trainable = False

        # Discriminators determines validity of translated images
#        valid_A,v_A2 = self.d_A(fake_A)
#        valid_B,v_B2 = self.d_B(fake_B)

        valid_A = self.d_A(fake_A)
        valid_B = self.d_B(fake_B)

        ssimLoss = DSSIMObjective()
        self.combined = Model(inputs=[img_A, img_B],
                              outputs=[ valid_A, valid_B,
                                        reconstr_A, reconstr_B,
                                        img_A_id, img_B_id ])
        self.combined.compile(loss=[lse, lse,
                                    ssimLoss, ssimLoss,
                                    'mae', 'mae'],
                            loss_weights=[  2, 2,
                                            self.lambda_cycle/5, self.lambda_cycle/5,
                                            0.1*self.lambda_id, 0.1*self.lambda_id ],
                            optimizer=optimizer)

    def build_generator(self):
        """U-Net Generator"""

        def conv2d(layer_input, filters, f_size=4):
            """Layers used during downsampling"""
            d = Conv2D(filters, kernel_size=f_size, strides=2, padding='same')(layer_input)
            d = LeakyReLU(alpha=0.2)(d)
            d = InstanceNormalization()(d)
            return d

        def deconv2d(layer_input, skip_input, filters, f_size=4, dropout_rate=0.1):
            """Layers used during upsampling"""
            u = UpSampling2D(size=2)(layer_input)
            u = Conv2D(filters, kernel_size=f_size, strides=1, padding='same', activation='relu')(u)
            if dropout_rate:
                u = Dropout(dropout_rate)(u)
            u = InstanceNormalization()(u)
            u = Concatenate()([u, skip_input])
            return u

        # Image input
        d0 = Input(shape=self.img_shape)

        # Downsampling
        d1 = conv2d(d0, self.gf)
        d2 = conv2d(d1, self.gf*2)
        d3 = conv2d(d2, self.gf*4)
        d4 = conv2d(d3, self.gf*8)

        # Upsampling•
        u1 = deconv2d(d4, d3, self.gf*4)
        u2 = deconv2d(u1, d2, self.gf*2)
        u3 = deconv2d(u2, d1, self.gf)

        u4 = UpSampling2D(size=2)(u3)
#        u4 = ReflectionPadding2D((3, 3))(u3)
        output_img = Conv2D(self.channels, kernel_size=4, strides=1, padding='same', activation='tanh')(u4)

        return Model(d0, output_img)

    def build_discriminator(self):
        def d_layer(layer_input, filters, f_size=4, normalization=True):
            """Discriminator layer"""
            d = Conv2D(filters, kernel_size=f_size, strides=2, padding='same')(layer_input)
            d = LeakyReLU(alpha=0.2)(d)
            if normalization:
                d = InstanceNormalization()(d)
            return d

        img = Input(shape=self.img_shape)

        d1 = d_layer(img, self.df, normalization=False)
        d2 = d_layer(d1, self.df*2)
        d3 = d_layer(d2, self.df*4)
        d4 = d_layer(d3, self.df*8)

        validity = Conv2D(1, kernel_size=4, strides=1, padding='same')(d4)
        return Model(img, validity)
    def train(self, epochs, epochStart=0,iteration=0, batch_size=1, sample_interval=50):
        start_time = datetime.datetime.now()

        # Adversarial loss ground truths
        valid = np.ones((batch_size,) + self.disc_patch)
        fake = np.zeros((batch_size,) + self.disc_patch)
        v2 = np.ones(batch_size)
        f2= np.zeros(batch_size)

        for epoch in range(epochStart,epochs):
            acc=[]

            for batch_i, (imgs_A, imgs_B) in enumerate(self.data_loader.load_batch(batch_size),iteration):

                # ----------------------
                #  Train Discriminators
                # ----------------------

                # Translate images to opposite domain
                fake_B = self.g_AB.predict(imgs_A)
                fake_A = self.g_BA.predict(imgs_B)

                # Train the discriminators (original images = real / translated = Fake)
                if( (batch_i%6)==0):
                    dA_loss_real = self.d_A.train_on_batch(imgs_A, [fake])
                    dA_loss_fake = self.d_A.train_on_batch(fake_A, [valid])
                    dB_loss_real = self.d_B.train_on_batch(imgs_B, [fake])
                    dB_loss_fake = self.d_B.train_on_batch(fake_B, [valid])
                else:
                    dA_loss_real = self.d_A.train_on_batch(imgs_A, [valid])
                    dA_loss_fake = self.d_A.train_on_batch(fake_A, [fake])
                    dB_loss_real = self.d_B.train_on_batch(imgs_B, [valid])
                    dB_loss_fake = self.d_B.train_on_batch(fake_B, [fake])
                dA_loss = 0.5 * np.add(dA_loss_real, dA_loss_fake)
                dB_loss = 0.5 * np.add(dB_loss_real, dB_loss_fake)
                
                # Total disciminator loss
                d_loss = 0.5 * np.add(dA_loss, dB_loss)
                # ------------------
                #  Train Generators
                # ------------------

                # Train the generators

                g_loss = self.combined.train_on_batch([imgs_A, imgs_B],
                                    [valid, valid,
                                    imgs_A, imgs_B,
                                    imgs_A, imgs_B])
                elapsed_time = datetime.datetime.now() - start_time

                print ("[Epoch %d/%d] [Batch %d/%d] [D loss: %f, acc: %3d%%] [G loss: %05f, adv: %05f, recon: %05f, id: %05f] time: %s " \
                                                                        % ( epoch, epochs,
                                                                            batch_i, self.data_loader.n_batches,
                                                                            d_loss[0], 100*d_loss[1],
                                                                            g_loss[0],
                                                                            np.mean(g_loss[1:3]),
                                                                            np.mean(g_loss[3:5]),
                                                                            np.mean(g_loss[5:6]),
                                                                            elapsed_time))

                acc.append(d_loss[1]*100)
                # If at save interval => save generated image samples
                if batch_i % sample_interval == 0:
                    self.sample_images(epoch, batch_i)


    def sample_images(self, epoch, batch_i):

        r, c = 6, 10

        imgs_A = self.data_loader.load_data(domain="A", dataName='TestCarlDatabaseFace', batch_size=c, is_testing=True)
        imgs_B = self.data_loader.load_data(domain="B", dataName='TestCarlDatabaseFace', batch_size=c, is_testing=True)

        # Translate images to the other domain

        fake_B = self.g_AB.predict(imgs_A)
        fake_A = self.g_BA.predict(imgs_B)

        # Translate back to original domain
        reconstr_A = self.g_BA.predict(fake_B)
        reconstr_B = self.g_AB.predict(fake_A)

        gen_imgs = np.concatenate([imgs_A, fake_B, reconstr_A, imgs_B, fake_A, reconstr_B])
        # Rescale images 0 - 1
        gen_imgs = 0.5 * gen_imgs + 0.5

        cnt = 0
        for i in range(r):
            for j in range(c):
                img = gen_imgs[cnt]
                if cnt==0 or (cnt%c)==0:
                    cImg = img
                else:
                    cImg = np.hstack((cImg, img))

                cnt += 1
            if i==0:
                vImg = cImg
            else:
                vImg = np.vstack((vImg,cImg))
        import cv2

        cv2.imshow('',cv2.resize(vImg,(640,384)))
        t = cv2.waitKey(1)
        ganType = 'carl'
        os.makedirs(ganType +'/img' , exist_ok=True)
        svStr = "%s/img/%03d_%d.png" % (ganType,epoch, batch_i)

        cv2.imwrite(svStr,vImg*255)
        os.makedirs(ganType +'/model', exist_ok=True)
        if (batch_i % 10000)==0:
            self.combined.save(ganType +"/model//model_"+str(epoch)+"-"+str(batch_i)+".hdf5")

        plt.close()


if __name__ == '__main__':
    import time
    tic = time.time()
    gan = CycleGAN()
    toc = time.time()
    es = '1'
    it='0'
    if es!='0' or it!='0':
        gan.combined.load_weights('carl\\model\\model_'+es+'-0.hdf5')
    print('Model Creation: ' + str(toc-tic))
    gan.train(epochs=1000, epochStart=int(es),iteration=int(it),batch_size=1, sample_interval=100)