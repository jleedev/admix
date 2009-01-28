class LocusFile < ActiveRecord::Base
  belongs_to :user
  has_many :marker
  validates_presence_of :name,:population
  validates_numericality_of :population
end
