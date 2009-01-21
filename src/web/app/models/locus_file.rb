class LocusFile < ActiveRecord::Base
  belongs_to :user
  has_many :marker
end
