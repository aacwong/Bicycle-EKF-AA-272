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

#import "FITDate.h"
#import "FITMessage.h"
#import "FITTypes.h"

NS_ASSUME_NONNULL_BEGIN

@interface FITVideoClipMesg : FITMessage
- (id)init;
// ClipNumber 
- (BOOL)isClipNumberValid;
- (FITUInt16)getClipNumber;
- (void)setClipNumber:(FITUInt16)clipNumber;
// StartTimestamp 
- (BOOL)isStartTimestampValid;
- (FITDate *)getStartTimestamp;
- (void)setStartTimestamp:(FITDate *)startTimestamp;
// StartTimestampMs 
- (BOOL)isStartTimestampMsValid;
- (FITUInt16)getStartTimestampMs;
- (void)setStartTimestampMs:(FITUInt16)startTimestampMs;
// EndTimestamp 
- (BOOL)isEndTimestampValid;
- (FITDate *)getEndTimestamp;
- (void)setEndTimestamp:(FITDate *)endTimestamp;
// EndTimestampMs 
- (BOOL)isEndTimestampMsValid;
- (FITUInt16)getEndTimestampMs;
- (void)setEndTimestampMs:(FITUInt16)endTimestampMs;
// ClipStart 
- (BOOL)isClipStartValid;
- (FITUInt32)getClipStart;
- (void)setClipStart:(FITUInt32)clipStart;
// ClipEnd 
- (BOOL)isClipEndValid;
- (FITUInt32)getClipEnd;
- (void)setClipEnd:(FITUInt32)clipEnd;

@end

NS_ASSUME_NONNULL_END
