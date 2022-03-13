////////////////////////////////////////////////////////////////////////////////
// The following FIT Protocol software provided may be used with FIT protocol
// devices only and remains the copyrighted property of Garmin Canada Inc.
// The software is being provided on an "as-is" basis and as an accommodation,
// and therefore all warranties, representations, or guarantees of any kind
// (whether express, implied or statutory) including, without limitation,
// warranties of merchantability, non-infringement, or fitness for a particular
// purpose, are specifically disclaimed.
//
// Copyright 2021 Garmin International, Inc.
////////////////////////////////////////////////////////////////////////////////
// ****WARNING****  This file is auto-generated!  Do NOT edit this file.
// Profile Version = 21.67Release
// Tag = production/akw/21.67.00-0-gd790f76b
////////////////////////////////////////////////////////////////////////////////


#import <Foundation/Foundation.h>

#import "FITMessage.h"
#import "FITTypes.h"

NS_ASSUME_NONNULL_BEGIN

@interface FITAntChannelIdMesg : FITMessage
- (id)init;
// ChannelNumber 
- (BOOL)isChannelNumberValid;
- (FITUInt8)getChannelNumber;
- (void)setChannelNumber:(FITUInt8)channelNumber;
// DeviceType 
- (BOOL)isDeviceTypeValid;
- (FITUInt8z)getDeviceType;
- (void)setDeviceType:(FITUInt8z)deviceType;
// DeviceNumber 
- (BOOL)isDeviceNumberValid;
- (FITUInt16z)getDeviceNumber;
- (void)setDeviceNumber:(FITUInt16z)deviceNumber;
// TransmissionType 
- (BOOL)isTransmissionTypeValid;
- (FITUInt8z)getTransmissionType;
- (void)setTransmissionType:(FITUInt8z)transmissionType;
// DeviceIndex 
- (BOOL)isDeviceIndexValid;
- (FITDeviceIndex)getDeviceIndex;
- (void)setDeviceIndex:(FITDeviceIndex)deviceIndex;

@end

NS_ASSUME_NONNULL_END
