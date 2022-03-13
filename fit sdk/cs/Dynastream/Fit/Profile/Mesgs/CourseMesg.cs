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
    /// Implements the Course profile message.
    /// </summary>
    public class CourseMesg : Mesg
    {
        #region Fields
        #endregion

        /// <summary>
        /// Field Numbers for <see cref="CourseMesg"/>
        /// </summary>
        public sealed class FieldDefNum
        {
            public const byte Sport = 4;
            public const byte Name = 5;
            public const byte Capabilities = 6;
            public const byte SubSport = 7;
            public const byte Invalid = Fit.FieldNumInvalid;
        }

        #region Constructors
        public CourseMesg() : base(Profile.GetMesg(MesgNum.Course))
        {
        }

        public CourseMesg(Mesg mesg) : base(mesg)
        {
        }
        #endregion // Constructors

        #region Methods
        ///<summary>
        /// Retrieves the Sport field</summary>
        /// <returns>Returns nullable Sport enum representing the Sport field</returns>
        public Sport? GetSport()
        {
            object obj = GetFieldValue(4, 0, Fit.SubfieldIndexMainField);
            Sport? value = obj == null ? (Sport?)null : (Sport)obj;
            return value;
        }

        /// <summary>
        /// Set Sport field</summary>
        /// <param name="sport_">Nullable field value to be set</param>
        public void SetSport(Sport? sport_)
        {
            SetFieldValue(4, 0, sport_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the Name field</summary>
        /// <returns>Returns byte[] representing the Name field</returns>
        public byte[] GetName()
        {
            byte[] data = (byte[])GetFieldValue(5, 0, Fit.SubfieldIndexMainField);
            return data.Take(data.Length - 1).ToArray();
        }

        ///<summary>
        /// Retrieves the Name field</summary>
        /// <returns>Returns String representing the Name field</returns>
        public String GetNameAsString()
        {
            byte[] data = (byte[])GetFieldValue(5, 0, Fit.SubfieldIndexMainField);
            return data != null ? Encoding.UTF8.GetString(data, 0, data.Length - 1) : null;
        }

        ///<summary>
        /// Set Name field</summary>
        /// <param name="name_"> field value to be set</param>
        public void SetName(String name_)
        {
            byte[] data = Encoding.UTF8.GetBytes(name_);
            byte[] zdata = new byte[data.Length + 1];
            data.CopyTo(zdata, 0);
            SetFieldValue(5, 0, zdata, Fit.SubfieldIndexMainField);
        }

        
        /// <summary>
        /// Set Name field</summary>
        /// <param name="name_">field value to be set</param>
        public void SetName(byte[] name_)
        {
            SetFieldValue(5, 0, name_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the Capabilities field</summary>
        /// <returns>Returns nullable uint representing the Capabilities field</returns>
        public uint? GetCapabilities()
        {
            Object val = GetFieldValue(6, 0, Fit.SubfieldIndexMainField);
            if(val == null)
            {
                return null;
            }

            return (Convert.ToUInt32(val));
            
        }

        /// <summary>
        /// Set Capabilities field</summary>
        /// <param name="capabilities_">Nullable field value to be set</param>
        public void SetCapabilities(uint? capabilities_)
        {
            SetFieldValue(6, 0, capabilities_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the SubSport field</summary>
        /// <returns>Returns nullable SubSport enum representing the SubSport field</returns>
        public SubSport? GetSubSport()
        {
            object obj = GetFieldValue(7, 0, Fit.SubfieldIndexMainField);
            SubSport? value = obj == null ? (SubSport?)null : (SubSport)obj;
            return value;
        }

        /// <summary>
        /// Set SubSport field</summary>
        /// <param name="subSport_">Nullable field value to be set</param>
        public void SetSubSport(SubSport? subSport_)
        {
            SetFieldValue(7, 0, subSport_, Fit.SubfieldIndexMainField);
        }
        
        #endregion // Methods
    } // Class
} // namespace
