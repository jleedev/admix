class CreateMarkers < ActiveRecord::Migration
  def self.up
    create_table :markers do |t|
      t.string :name, :limit => 16
      t.belongs_to :locus_file

      t.timestamps
    end
  end

  def self.down
    drop_table :markers
  end
end
