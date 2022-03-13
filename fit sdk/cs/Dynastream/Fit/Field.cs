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

using System.Collections.Generic;
using System.Linq;

namespace Dynastream.Fit
{
    public class Field
        : FieldBase
    {
        #region Fields
        private string name;
        private byte type;
        private double scale;
        private double offset;
        private string units;
        private bool isAccumulated;
        private Profile.Type profileType;

        internal List<Subfield> subfields = new List<Subfield>();
        internal List<FieldComponent> components = new List<FieldComponent>();
        #endregion

        #region Properties
        public override string Name
        {
            get
            {
                return name;
            }
        }

        public byte Num { get; set; }

        public override byte Type
        {
            get
            {
                return type;
            }
        }

        public override double Scale
        {
            get
            {
                return scale;
            }
        }

        public override double Offset
        {
            get
            {
                return offset;
            }
        }

        public override string Units
        {
            get
            {
                return units;
            }
        }

        public bool IsAccumulated
        {
            get
            {
                return isAccumulated;
            }
        }

        public Profile.Type ProfileType
        {
            get
            {
                return profileType;
            }
        }

        public bool IsExpandedField { get; set; }
        #endregion

        #region Constructors
        public Field(Field other)
            : base(other)
        {
            if (other == null)
            {
                this.name = "unknown";
                this.Num = Fit.FieldNumInvalid;
                this.type = 0;
                this.scale = 1f;
                this.offset = 0f;
                this.units = "";
                this.isAccumulated = false;
                this.profileType = Profile.Type.Enum;
                this.IsExpandedField = false;
                return;
            }

            this.name = other.Name;
            this.Num = other.Num;
            this.type = other.Type;
            this.scale = other.Scale;
            this.offset = other.Offset;
            this.units = other.units;
            this.isAccumulated = other.isAccumulated;
            this.profileType = other.profileType;
            this.IsExpandedField = other.IsExpandedField;

            foreach (Subfield subfield in other.subfields)
            {
                this.subfields.Add(new Subfield(subfield));
            }
            foreach (FieldComponent component in other.components)
            {
                this.components.Add(new FieldComponent(component));
            }
        }

        internal Field(string name, byte num, byte type, double scale, double offset, string units, bool accumulated, Profile.Type profileType)
        {
            this.name = name;
            this.Num = num;
            this.type = type;
            this.scale = scale;
            this.offset = offset;
            this.units = units;
            this.isAccumulated = accumulated;
            this.profileType = profileType;
            this.IsExpandedField = false;
        }

        internal Field(byte num, byte type)
            : this("unknown", num, type, 1.0d, 0.0d, "", false, Profile.Type.NumTypes)
        {
        }
        #endregion

        #region Methods

        internal void SetType(byte value)
        {
            type = value;
        }

        internal override Subfield GetSubfield(string subfieldName)
        {
            return subfields.FirstOrDefault(subfield => subfield.Name == subfieldName);
        }

        internal override Subfield GetSubfield(int subfieldIndex)
        {
            // SubfieldIndexActiveSubfield and SubfieldIndexMainField
            // will be out of this range
            if (subfieldIndex >= 0 && subfieldIndex < subfields.Count)
            {
                return subfields[subfieldIndex];
            }

            return null;

        }
        #endregion
    }
} // namespace
