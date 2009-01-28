class AddNameToLocusFile < ActiveRecord::Migration
  def self.up
    add_column :locus_files, :name, :string
  end

  def self.down
    remove_column :locus_files, :name
  end
end
