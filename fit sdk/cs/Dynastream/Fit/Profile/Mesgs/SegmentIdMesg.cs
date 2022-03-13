#region Copyright
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

#endregion

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using System.IO;
using System.Linq;

namespace Dynastream.Fit
{
    /// <summary>
    /// Implements the SegmentId profile message.
    /// </summary>
    public class SegmentIdMesg : Mesg
    {
        #region Fields
        #endregion

        /// <summary>
        /// Field Numbers for <see cref="SegmentIdMesg"/>
        /// </summary>
        public sealed class FieldDefNum
        {
            public const byte Name = 0;
            public const byte Uuid = 1;
            public const byte Sport = 2;
            public const byte Enabled = 3;
            public const byte UserProfilePrimaryKey = 4;
            public const byte DeviceId = 5;
            public const byte DefaultRaceLeader = 6;
            public const byte DeleteStatus = 7;
            public const byte SelectionType = 8;
            public const byte Invalid = Fit.FieldNumInvalid;
        }

        #region Constructors
        public SegmentIdMesg() : base(Profile.GetMesg(MesgNum.SegmentId))
        {
        }

        public SegmentIdMesg(Mesg mesg) : base(mesg)
        {
        }
        #endregion // Constructors

        #region Methods
        ///<summary>
        /// Retrieves the Name field
        /// Comment: Friendly name assigned to segment</summary>
        /// <returns>Returns byte[] representing the Name field</returns>
        public byte[] GetName()
        {
            byte[] data = (byte[])GetFieldValue(0, 0, Fit.SubfieldIndexMainField);
            return data.Take(data.Length - 1).ToArray();
        }

        ///<summary>
        /// Retrieves the Name field
        /// Comment: Friendly name assigned to segment</summary>
        /// <returns>Returns String representing the Name field</returns>
        public String GetNameAsString()
        {
            byte[] data = (byte[])GetFieldValue(0, 0, Fit.SubfieldIndexMainField);
            return data != null ? Encoding.UTF8.GetString(data, 0, data.Length - 1) : null;
        }

        ///<summary>
        /// Set Name field
        /// Comment: Friendly name assigned to segment</summary>
        /// <param name="name_"> field value to be set</param>
        public void SetName(String name_)
        {
            byte[] data = Encoding.UTF8.GetBytes(name_);
            byte[] zdata = new byte[data.Length + 1];
            data.CopyTo(zdata, 0);
            SetFieldValue(0, 0, zdata, Fit.SubfieldIndexMainField);
        }

        
        /// <summary>
        /// Set Name field
        /// Comment: Friendly name assigned to segment</summary>
        /// <param name="name_">field value to be set</param>
        public void SetName(byte[] name_)
        {
            SetFieldValue(0, 0, name_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the Uuid field
        /// Comment: UUID of the segment</summary>
        /// <returns>Returns byte[] representing the Uuid field</returns>
        public byte[] GetUuid()
        {
            byte[] data = (byte[])GetFieldValue(1, 0, Fit.SubfieldIndexMainField);
            return data.Take(data.Length - 1).ToArray();
        }

        ///<summary>
        /// Retrieves the Uuid field
        /// Comment: UUID of the segment</summary>
        /// <returns>Returns String representing the Uuid field</returns>
        public String GetUuidAsString()
        {
            byte[] data = (byte[])GetFieldValue(1, 0, Fit.SubfieldIndexMainField);
            return data != null ? Encoding.UTF8.GetString(data, 0, data.Length - 1) : null;
        }

        ///<summary>
        /// Set Uuid field
        /// Comment: UUID of the segment</summary>
        /// <param name="uuid_"> field value to be set</param>
        public void SetUuid(String uuid_)
        {
            byte[] data = Encoding.UTF8.GetBytes(uuid_);
            byte[] zdata = new byte[data.Length + 1];
            data.CopyTo(zdata, 0);
            SetFieldValue(1, 0, zdata, Fit.SubfieldIndexMainField);
        }

        
        /// <summary>
        /// Set Uuid field
        /// Comment: UUID of the segment</summary>
        /// <param name="uuid_">field value to be set</param>
        public void SetUuid(byte[] uuid_)
        {
            SetFieldValue(1, 0, uuid_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the Sport field
        /// Comment: Sport associated with the segment</summary>
        /// <returns>Returns nullable Sport enum representing the Sport field</returns>
        public Sport? GetSport()
        {
            object obj = GetFieldValue(2, 0, Fit.SubfieldIndexMainField);
            Sport? value = obj == null ? (Sport?)null : (Sport)obj;
            return value;
        }

        /// <summary>
        /// Set Sport field
        /// Comment: Sport associated with the segment</summary>
        /// <param name="sport_">Nullable field value to be set</param>
        public void SetSport(Sport? sport_)
        {
            SetFieldValue(2, 0, sport_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the Enabled field
        /// Comment: Segment enabled for evaluation</summary>
        /// <returns>Returns nullable Bool enum representing the Enabled field</returns>
        public Bool? GetEnabled()
        {
            object obj = GetFieldValue(3, 0, Fit.SubfieldIndexMainField);
            Bool? value = obj == null ? (Bool?)null : (Bool)obj;
            return value;
        }

        /// <summary>
        /// Set Enabled field
        /// Comment: Segment enabled for evaluation</summary>
        /// <param name="enabled_">Nullable field value to be set</param>
        public void SetEnabled(Bool? enabled_)
        {
            SetFieldValue(3, 0, enabled_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the UserProfilePrimaryKey field
        /// Comment: Primary key of the user that created the segment</summary>
        /// <returns>Returns nullable uint representing the UserProfilePrimaryKey field</returns>
        public uint? GetUserProfilePrimaryKey()
        {
            Object val = GetFieldValue(4, 0, Fit.SubfieldIndexMainField);
            if(val == null)
            {
                return null;
            }

            return (Convert.ToUInt32(val));
            
        }

        /// <summary>
        /// Set UserProfilePrimaryKey field
        /// Comment: Primary key of the user that created the segment</summary>
        /// <param name="userProfilePrimaryKey_">Nullable field value to be set</param>
        public void SetUserProfilePrimaryKey(uint? userProfilePrimaryKey_)
        {
            SetFieldValue(4, 0, userProfilePrimaryKey_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the DeviceId field
        /// Comment: ID of the device that created the segment</summary>
        /// <returns>Returns nullable uint representing the DeviceId field</returns>
        public uint? GetDeviceId()
        {
            Object val = GetFieldValue(5, 0, Fit.SubfieldIndexMainField);
            if(val == null)
            {
                return null;
            }

            return (Convert.ToUInt32(val));
            
        }

        /// <summary>
        /// Set DeviceId field
        /// Comment: ID of the device that created the segment</summary>
        /// <param name="deviceId_">Nullable field value to be set</param>
        public void SetDeviceId(uint? deviceId_)
        {
            SetFieldValue(5, 0, deviceId_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the DefaultRaceLeader field
        /// Comment: Index for the Leader Board entry selected as the default race participant</summary>
        /// <returns>Returns nullable byte representing the DefaultRaceLeader field</returns>
        public byte? GetDefaultRaceLeader()
        {
            Object val = GetFieldValue(6, 0, Fit.SubfieldIndexMainField);
            if(val == null)
            {
                return null;
            }

            return (Convert.ToByte(val));
            
        }

        /// <summary>
        /// Set DefaultRaceLeader field
        /// Comment: Index for the Leader Board entry selected as the default race participant</summary>
        /// <param name="defaultRaceLeader_">Nullable field value to be set</param>
        public void SetDefaultRaceLeader(byte? defaultRaceLeader_)
        {
            SetFieldValue(6, 0, defaultRaceLeader_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the DeleteStatus field
        /// Comment: Indicates if any segments should be deleted</summary>
        /// <returns>Returns nullable SegmentDeleteStatus enum representing the DeleteStatus field</returns>
        public SegmentDeleteStatus? GetDeleteStatus()
        {
            object obj = GetFieldValue(7, 0, Fit.SubfieldIndexMainField);
            SegmentDeleteStatus? value = obj == null ? (SegmentDeleteStatus?)null : (SegmentDeleteStatus)obj;
            return value;
        }

        /// <summary>
        /// Set DeleteStatus field
        /// Comment: Indicates if any segments should be deleted</summary>
        /// <param name="deleteStatus_">Nullable field value to be set</param>
        public void SetDeleteStatus(SegmentDeleteStatus? deleteStatus_)
        {
            SetFieldValue(7, 0, deleteStatus_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the SelectionType field
        /// Comment: Indicates how the segment was selected to be sent to the device</summary>
        /// <returns>Returns nullable SegmentSelectionType enum representing the SelectionType field</returns>
        public SegmentSelectionType? GetSelectionType()
        {
            object obj = GetFieldValue(8, 0, Fit.SubfieldIndexMainField);
            SegmentSelectionType? value = obj == null ? (SegmentSelectionType?)null : (SegmentSelectionType)obj;
            return value;
        }

        /// <summary>
        /// Set SelectionType field
        /// Comment: Indicates how the segment was selected to be sent to the device</summary>
        /// <param name="selectionType_">Nullable field value to be set</param>
        public void SetSelectionType(SegmentSelectionType? selectionType_)
        {
            SetFieldValue(8, 0, selectionType_, Fit.SubfieldIndexMainField);
        }
        
        #endregion // Methods
    } // Class
} // namespace
