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

#import <Foundation/Foundation.h>

#import "FITMessage.h"
#import "FITDeveloperFieldDefinition.h"

NS_ASSUME_NONNULL_BEGIN

@class FITDecoder;
@protocol FITMesgDelegate <NSObject>
- (void)onMesg:(FITMessage*)mesg;
@end

@protocol FITMesgDefinitionDelegate  <NSObject>
- (void)onMesgDefinition:(FITMessage*)mesg;
@end

@protocol FITDeveloperFieldDefinitionDelegate  <NSObject>
- (void)onDeveloperFieldDefinition:(FITDeveloperFieldDefinition*)definition;
@end

@interface FITDecoder : NSObject
@property (strong) id<FITMesgDelegate>mesgDelegate;
@property (strong) id<FITMesgDefinitionDelegate>mesgDefinitionDelegate;
@property (strong) id<FITDeveloperFieldDefinitionDelegate>developerFieldDefinitionDelegate;
- (instancetype)init;
- (BOOL)isFIT:(NSString*)filename;
- (BOOL)checkIntegrity:(NSString*)filename;
- (BOOL)decodeFile:(NSString *)filename;

@end

NS_ASSUME_NONNULL_END
